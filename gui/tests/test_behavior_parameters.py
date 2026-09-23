"""Parameter discovery, validation and temporary-sketch integration tests."""
from dataclasses import asdict
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

from gui import behavior as gui
from gui.parameters import apply_overrides, parameters_in, prepare_sketch, source_hash, validate_value


SOURCE = b'''// User settings
boolean ENL = true; // enforce no licking
unsigned long UnitRewardSize = 24; // calibration in ms
unsigned long SmallRewardSize = 2 * UnitRewardSize;
unsigned long BigRewardSize = 10 * UnitRewardSize;
// Outcome probability params
int PairProbRange[2] = {-1, -40};
float gain = 0.5;
// unsigned long Commented = 123;
/* unsigned long AlsoCommented = 123; */
// Params Initializtion
int TrialNum = 0;
void setup() { unsigned long UnitRewardSize = 99; }
void loop() {}
'''


class ParameterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.sketch = self.base / 'Reward'
        self.sketch.mkdir()
        self.main = self.sketch / 'Reward.ino'
        self.main.write_bytes(SOURCE)

    def test_only_supported_settings_not_runtime_locals_or_comments(self):
        parameters = parameters_in(SOURCE.decode())
        self.assertEqual([p.name for p in parameters],
                         ['ENL', 'UnitRewardSize', 'SmallRewardSize', 'BigRewardSize', 'PairProbRange', 'gain'])
        self.assertEqual(parameters[1].comment, 'calibration in ms')
        self.assertEqual(parameters[-2].group, 'Outcome probability params')
        self.assertEqual(parameters[-2].size, 2)

    def test_explicit_markers_are_supported_and_unmarked_sources_are_untouched(self):
        source = '// GUI parameters begin\nint count = 1;\n// GUI parameters end\nint state = 0;'
        self.assertEqual([p.name for p in parameters_in(source)], ['count'])
        self.assertEqual(parameters_in('int count = 1;'), [])
        self.assertEqual(parameters_in('// User settings\nint count = 1;'), [])

    def test_excludes_multideclarations_functions_conditionals_and_strings(self):
        source = '''// User settings
int a = 1, b = 2;
int value = random(0, 20);
int string = "// not a parameter";
#if SOMETHING
int conditional = 10;
#else
int conditional = 20;
#endif
void helper() {
int local = 2;
}
int duplicate = 1;
int duplicate = 2;
int actual = 3;
// Params Initialization
'''
        self.assertEqual([p.name for p in parameters_in(source)], ['actual'])

    def test_exact_span_patch_preserves_formulas_comments_and_crlf(self):
        source = SOURCE.replace(b'\n', b'\r\n')
        result = apply_overrides(source, {'UnitRewardSize': '30'}, source_hash(source))
        self.assertEqual(result, source.replace(b'UnitRewardSize = 24;', b'UnitRewardSize = 30;'))
        self.assertIn(b'2 * UnitRewardSize', result)
        self.assertIn(b'UnitRewardSize = 99;', result)

    def test_bool_array_and_formula_edits(self):
        edits = {'ENL': 'false', 'PairProbRange': '{1, 40}', 'BigRewardSize': '8 * UnitRewardSize'}
        result = apply_overrides(SOURCE, edits, source_hash(SOURCE)).decode()
        self.assertIn('boolean ENL = false;', result)
        self.assertIn('int PairProbRange[2] = {1, 40};', result)
        self.assertIn('BigRewardSize = 8 * UnitRewardSize;', result)

    def test_unit_requires_positive_integer_and_disallows_code(self):
        p = next(p for p in parameters_in(SOURCE.decode()) if p.name == 'UnitRewardSize')
        for value in ('', '0', '-1', '1.5', '1/2', '5; exit()', '4294967296', 'abc'):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'UnitRewardSize'):
                validate_value(p, value)

    def test_rejects_invalid_array_bool_and_expression_syntax(self):
        for name, value in [('ENL', 'yes'), ('PairProbRange', '{1}'), ('PairProbRange', '{1, 2, 3}'),
                            ('SmallRewardSize', '1; void loop() {}'), ('SmallRewardSize', 'random(2,3)'),
                            ('SmallRewardSize', '1 // 2'), ('gain', '__import__("os")'),
                            ('gain', '2 ** 3'), ('SmallRewardSize', '-5')]:
            with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                apply_overrides(SOURCE, {name: value}, source_hash(SOURCE))

    def test_cpp_numeric_suffixes_and_expressions_remain_intact(self):
        for value in ('2UL * UnitRewardSize', '(UnitRewardSize + 2) / 2', '0xFF', '5U'):
            with self.subTest(value=value):
                result = apply_overrides(SOURCE, {'SmallRewardSize': value}, source_hash(SOURCE))
                self.assertIn(value.encode(), result)

    def test_unknown_names_and_stale_edits_fail(self):
        with self.assertRaisesRegex(ValueError, 'not editable'):
            apply_overrides(SOURCE, {'TrialNum': '5'}, source_hash(SOURCE))
        with self.assertRaisesRegex(ValueError, 'source has changed'):
            apply_overrides(SOURCE + b'// edit\n', {'UnitRewardSize': '30'}, source_hash(SOURCE))

    def test_whole_sketch_is_copied_and_only_main_initializer_changes(self):
        (self.sketch / 'src').mkdir()
        (self.sketch / 'src' / 'helper.h').write_bytes(b'#define X 1\r\n')
        (self.sketch / 'Another.ino').write_bytes(b'void helper() {}')
        before = {p.relative_to(self.sketch): p.read_bytes() for p in self.sketch.rglob('*') if p.is_file()}
        staged = prepare_sketch(self.sketch, self.base / 'staged', {'UnitRewardSize': '30'}, source_hash(SOURCE))
        self.assertEqual(staged.name, self.sketch.name)
        self.assertEqual((staged / 'Reward.ino').read_bytes(), SOURCE.replace(b'= 24;', b'= 30;'))
        for relative, data in before.items():
            self.assertEqual((self.sketch / relative).read_bytes(), data)
            if relative.name != 'Reward.ino':
                self.assertEqual((staged / relative).read_bytes(), data)

    def test_upload_uses_edited_copy_and_preserves_original_on_compile_failure(self):
        service = gui.BehaviorService()
        settings = gui.Settings(sketch=str(self.main), reconnect=False,
                                parameter_overrides={'UnitRewardSize': '32'}, parameter_source_hash=source_hash(SOURCE))
        seen = []
        def compile_fails(args, timeout):
            staged = Path(args[-1])
            seen.append(staged)
            self.assertIn(b'UnitRewardSize = 32;', (staged / 'Reward.ino').read_bytes())
            raise RuntimeError('compile failed')
        with patch.object(gui, 'cli_executable', return_value='cli'), patch.object(service, '_run_cli', side_effect=compile_fails) as run:
            with self.assertRaisesRegex(RuntimeError, 'compile failed'):
                service.do_upload(settings)
        run.assert_called_once()
        self.assertEqual(self.main.read_bytes(), SOURCE)
        self.assertFalse(seen[0].exists())
        self.assertFalse(service.busy)

    def test_stale_upload_does_not_disconnect_or_compile(self):
        service = gui.BehaviorService()
        serial = Mock()
        service.ser = serial
        settings = gui.Settings(sketch=str(self.main), parameter_overrides={'UnitRewardSize': '32'},
                                parameter_source_hash='old-source')
        with patch.object(gui, 'cli_executable', return_value='cli'), patch.object(service, '_run_cli') as run:
            with self.assertRaisesRegex(ValueError, 'source has changed'):
                service.do_upload(settings)
        serial.close.assert_not_called()
        run.assert_not_called()

    def test_profile_roundtrip_and_old_profile_defaults(self):
        settings = gui.Settings(sketch=str(self.main), parameter_overrides={'UnitRewardSize': '30'},
                                parameter_source_hash=source_hash(SOURCE))
        restored = gui.Settings(**json.loads(json.dumps(asdict(settings))))
        self.assertEqual(restored, settings)
        self.assertIn(b'UnitRewardSize = 30;', apply_overrides(SOURCE, restored.parameter_overrides, restored.parameter_source_hash))
        old = gui.Settings(**{'port': 'COM5', 'baud': '9600'})
        self.assertEqual(old.parameter_overrides, {})


if __name__ == '__main__':
    unittest.main()
