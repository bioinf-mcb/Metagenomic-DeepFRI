import time


class ElapsedTimeHandler(object):
    def __init__(self, path):
        self.path = path
        self.start = time.time()

    def log(self, name):
        text = f"{name},{time.time() - self.start}\n"
        print(text[:-1])    # omitting newline
        with open(self.path, "a") as f:
            f.write(text)
        self.reset()

    def reset(self):
        self.start = time.time()
