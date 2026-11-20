import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

from gui.reward_receiver import RewardReceiver
from gui.behavior import BehaviorService

CONFIG = dict(name='Test', threshold=10, reset_threshold=8, hold_ms=100, cooldown_s=1,
              reward='small', max_rewards=2)


class ReceiverTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.send = Mock()
        self.ready = Mock(return_value=True)
        with patch('gui.reward_receiver.socket.socket'):
            self.receiver = RewardReceiver(5006, '127.0.0.1', self.temp.name, self.send, Mock(), self.ready)
        self.receiver.client = Mock()
        self.addCleanup(self.receiver.close)
        self.clock = patch('gui.reward_receiver.time.monotonic', return_value=10)
        self.time = self.clock.start()
        self.addCleanup(self.clock.stop)
        self.receiver._handle(dict(type='start', config=dict(CONFIG)))

    def reward(self, **kwargs):
        message = dict(type='reward', token=self.receiver.token, index=self.receiver.count + 1, value=11, sample_s=1)
        message.update(kwargs)
        self.receiver._handle(message)

    def test_duplicate_token_never_sends_twice(self):
        token = self.receiver.token
        self.reward()
        self.time.return_value = 11
        with self.assertRaises(ValueError):
            self.reward(token=token, index=1)
        self.send.assert_called_once_with('w')

    def test_expired_request_is_rejected(self):
        self.time.return_value = 11.01
        with self.assertRaises(ValueError):
            self.reward()
        self.send.assert_not_called()

    def test_cooldown_and_limit_enforced_at_receiver(self):
        self.reward()
        with self.assertRaises(ValueError):
            self.reward()
        self.assertEqual(self.send.call_count, 1)
        # Rejections are fatal over the network; directly restore a token for this unit check.
        self.receiver._reply()
        self.time.return_value = 11
        self.reward()
        self.time.return_value = 12
        with self.assertRaises(ValueError):
            self.reward()
        self.assertEqual(self.send.call_count, 2)

    def test_not_ready_and_invalid_signal_never_send(self):
        self.ready.return_value = False
        with self.assertRaises(ValueError):
            self.reward()
        self.ready.return_value = True
        with self.assertRaises(ValueError):
            self.reward(value=float('nan'))
        self.send.assert_not_called()

    def test_log_failure_prevents_serial_write(self):
        with patch.object(self.receiver, 'record', side_effect=OSError('disk full')):
            with self.assertRaises(OSError):
                self.reward()
        self.send.assert_not_called()

    def test_serial_write_failure_is_not_reported_as_sent(self):
        self.send.side_effect = OSError('USB disconnected')
        with self.assertRaises(OSError):
            self.reward()
        self.send.assert_called_once()
        self.assertEqual(self.receiver.count, 0)
        rows = [json.loads(line) for line in next(Path(self.temp.name).glob('*.jsonl')).read_text().splitlines()]
        self.assertEqual([row['event'] for row in rows], ['task_start', 'reward_request'])

    def test_raw_arduino_responses_are_persisted_after_task_end(self):
        self.receiver._drop_client('Stopped')
        self.receiver.serial_line('ACK SMALL_REWARD')
        rows = [json.loads(line) for line in next(Path(self.temp.name).glob('*.jsonl')).read_text().splitlines()]
        self.assertEqual(rows[-1]['line'], 'ACK SMALL_REWARD')

    def test_bad_message_drops_connection_without_reward(self):
        self.receiver.client.recv.return_value = b'{invalid}\n'
        self.receiver.poll()
        self.assertIsNone(self.receiver.client)
        self.send.assert_not_called()


class ServiceTests(unittest.TestCase):
    def test_disconnect_disables_receiver(self):
        service = BehaviorService()
        receiver = service.reward_receiver = Mock()
        service.do_disconnect()
        receiver.close.assert_called_once()
        self.assertIsNone(service.reward_receiver)

    def test_typed_stop_or_calibration_disables_receiver_before_serial_write(self):
        for command in ('x', 'c'):
            service = BehaviorService()
            order = []
            receiver = service.reward_receiver = Mock()
            receiver.close.side_effect = lambda: order.append('closed')
            service.ser = Mock()
            service.ser.write.side_effect = lambda data: (order.append(data), len(data))[1]
            service.do_send(command, 'None')
            self.assertEqual(order, ['closed', command.encode()])

    def test_legacy_sketch_needs_no_capability_banner_for_rewards(self):
        service = BehaviorService()
        service.ser = Mock()
        service.ser.write.side_effect = len
        service.do_send('w', 'None')
        service.ser.write.assert_called_once_with(b'w')
        self.assertFalse(service.manual_ready)


if __name__ == '__main__':
    unittest.main()
