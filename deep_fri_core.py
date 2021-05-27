# it feels to me like all models should be pre_initialized as runtime services
# so that for example no GPU allocation would take place during main_pipeline

class DeepFRI:
    def __init__(self):
        print("Init DeepFRI")

    def __call__(self, sequence, contact_map):
        print("DeepFRI CALL")


class DeepCNN:
    def __init__(self):
        print("Init DeepCNN")

    def __call__(self, *args, **kwargs):
        print("DeepCNN CALL")


class DeepLM:
    def __init__(self):
        print("Init DeepLM")

    def __call__(self, *args, **kwargs):
        print("DeepLM CALL")
        # s = <sequence of length = L>
	# z := LM(s) = (L, d)
	# z @ z.T (L, L) 
        return 0
