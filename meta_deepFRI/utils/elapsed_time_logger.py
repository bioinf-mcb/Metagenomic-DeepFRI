import time


class ElapsedTimeLogger(object):

    def __init__(self, path=None):
        self.path = path
        self.current_start = time.time()
        self.start = self.current_start

    def log(self, name):
        text = f"{name},{(time.time() - self.current_start):.10f}\n"
        self._append_to_file(text)
        self.reset()

    def reset(self):
        self.current_start = time.time()

    def log_total_time(self):
        text = f"total_time,{(time.time() - self.start):.10f}\n"
        self._append_to_file(text)

    def _append_to_file(self, text):
        print(text[:-1])  # omitting newline
        if self.path is not None:
            with open(self.path, "a") as f:
                f.write(text)
