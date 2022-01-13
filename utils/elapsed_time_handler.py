import time


class ElapsedTimeHandler(object):
    def __init__(self, path):
        self.path = path
        self.start = time.time()

    def log(self, name):
        with open(self.path, "a") as f:
            f.write(f"{name},{time.time() - self.start}\n")
        self.reset()

    def reset(self):
        self.start = time.time()
