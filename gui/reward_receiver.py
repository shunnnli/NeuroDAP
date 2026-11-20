"""Opt-in BrainClamp reward receiver. Poll only on the Arduino's service thread."""
import json
import math
from pathlib import Path
import secrets
import socket
import time
import uuid


class RewardReceiver:
    def __init__(self, port, allowed_host, log_dir, send_reward, log, ready, bind_host='0.0.0.0'):
        self.allowed_host = allowed_host
        self.log_dir = Path(log_dir)
        self.send_reward, self.log, self.ready = send_reward, log, ready
        self.client = self.journal = None
        self.buffer = bytearray()
        self.config = None
        self.token = None
        self.count = 0
        self.listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            if hasattr(socket, 'SO_EXCLUSIVEADDRUSE'):
                self.listener.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
            else:
                self.listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            self.listener.bind((bind_host, port))
            self.listener.listen(1)
            self.listener.setblocking(False)
        except Exception:
            self.listener.close()
            raise
        self.port = self.listener.getsockname()[1]

    def record(self, event, **data):
        if self.journal:
            self.journal.write(json.dumps(dict(time_unix_s=time.time(), event=event, **data), allow_nan=False) + '\n')
            self.journal.flush()

    def serial_line(self, line):
        self.record('arduino_rx', line=line)

    def _drop_client(self, reason):
        client, self.client = self.client, None
        self.config = None
        self.token = None
        self.buffer.clear()
        if client:
            client.close()
            self.record('task_end', reason=reason, reward_commands_sent=self.count)
            self.log(f'[reward task ended] {reason}; {self.count} reward commands sent')

    def close(self):
        try:
            self._drop_client('Receiver disabled / Arduino disconnected')
        finally:
            self.listener.close()
            if self.journal:
                self.journal.close()
                self.journal = None

    def _reply(self, **data):
        self.token = secrets.token_hex(16)
        self.token_time = time.monotonic()
        self.client.sendall((json.dumps(dict(status='ok', token=self.token, **data)) + '\n').encode())

    @staticmethod
    def validate_config(config):
        if not isinstance(config, dict) or not isinstance(config.get('name'), str) or not 1 <= len(config['name']) <= 100:
            raise ValueError('Invalid task name.')
        for key in ('threshold', 'reset_threshold', 'hold_ms', 'cooldown_s'):
            value = config.get(key)
            if type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError(f'Invalid {key}.')
        if config['reset_threshold'] >= config['threshold'] or not 0 <= config['hold_ms'] <= 60000:
            raise ValueError('Invalid thresholds or hold time.')
        if not 0.1 <= config['cooldown_s'] <= 86400:
            raise ValueError('Invalid cooldown.')
        if type(config.get('max_rewards')) is not int or not 1 <= config['max_rewards'] <= 100000:
            raise ValueError('Invalid reward limit.')
        if config.get('reward') not in ('small', 'large'):
            raise ValueError('Invalid reward size.')
        return config

    def _handle(self, message):
        if not isinstance(message, dict):
            raise ValueError('Expected a task message.')
        kind = message.get('type')
        now = time.monotonic()
        if not self.ready():
            raise ValueError('Arduino is unavailable or busy.')
        if kind == 'start':
            if self.config is not None:
                raise ValueError('A task is already running on this connection.')
            config = self.validate_config(message.get('config'))
            if self.journal:
                self.journal.close()
                self.journal = None
            self.log_dir.mkdir(parents=True, exist_ok=True)
            path = self.log_dir / (time.strftime('%Y%m%d-%H%M%S-') + uuid.uuid4().hex + '.jsonl')
            self.journal = path.open('x', encoding='utf-8')
            self.count = 0
            self.last_reward = -math.inf
            self.record('task_start', config=config, peer=self.allowed_host)
            self.config = config
            self.log(f'[reward task started] {config["name"]}; log: {path}')
            self._reply()
            return
        if self.config is None or message.get('token') != self.token or now - self.token_time > 1:
            raise ValueError('Expired, duplicate, or inactive task request.')
        # Consume before any side effect. A command is never retried.
        self.token = None
        if kind == 'heartbeat':
            self._reply()
        elif kind == 'end':
            self._reply()
            self._drop_client(str(message.get('reason', 'Task stopped'))[:200])
        elif kind == 'reward':
            c = self.config
            if type(message.get('index')) is not int or message['index'] != self.count + 1:
                raise ValueError('Duplicate or out-of-order reward.')
            for key in ('value', 'sample_s'):
                if type(message.get(key)) not in (int, float) or not math.isfinite(message[key]):
                    raise ValueError('Invalid photometry sample.')
            if message['value'] <= c['threshold']:
                raise ValueError('Signal is not above the task threshold.')
            if self.count >= c['max_rewards'] or now - self.last_reward < c['cooldown_s']:
                raise ValueError('Reward limit or cooldown reached.')
            command = 'w' if c['reward'] == 'small' else 'W'
            self.record('reward_request', index=message['index'], value=message['value'],
                        sample_s=message['sample_s'], command=command)
            self.send_reward(command)
            self.count += 1
            self.last_reward = time.monotonic()
            self.record('reward_serial_sent', index=self.count, command=command)
            self.log(f'[reward task] #{self.count}: DA={message["value"]:.6g}, sent {command!r}')
            # This acknowledges a serial write, not physical water delivery.
            self._reply(index=self.count)
        else:
            raise ValueError('Unknown task request.')

    def poll(self):
        if self.client is None:
            try:
                client, peer = self.listener.accept()
            except BlockingIOError:
                return
            if peer[0] != self.allowed_host:
                client.close()
                return
            client.setblocking(False)
            self.client = client
            self.last_received = time.monotonic()
            self.buffer.clear()
        try:
            try:
                packet = self.client.recv(16384)
            except BlockingIOError:
                packet = None
            if packet == b'':
                self._drop_client('BrainClamp stopped or disconnected')
                return
            if packet:
                self.buffer.extend(packet)
                if len(self.buffer) > 16384:
                    raise ValueError('Task message too large.')
            if b'\n' in self.buffer:
                line, _, rest = self.buffer.partition(b'\n')
                self.buffer = bytearray(rest)
                self._handle(json.loads(line))
                self.last_received = time.monotonic()
            if self.client and time.monotonic() - self.last_received > 1.5:
                raise ValueError('BrainClamp heartbeat timed out.')
        except Exception as exc:
            self.record('task_error', error=str(exc))
            if self.client:
                try:
                    self.client.sendall((json.dumps(dict(status='error', error=str(exc))) + '\n').encode())
                except OSError:
                    pass
            self._drop_client(str(exc))
