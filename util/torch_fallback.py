import sys
import builtins  # <-- We use this to make torch completely global

try:
    import torch
    import torch.nn.functional as F
except ImportError:

    class MetaMock(type):
        def __getattr__(cls, name):
            return cls

        def __call__(cls, *args, **kwargs):
            return cls

    class MockTorch(metaclass=MetaMock):
        pass

    mock_instance = MockTorch

    # 1. Fix for files that DO execute "import torch"
    sys.modules["torch"] = mock_instance
    sys.modules["torch.nn.functional"] = mock_instance

    # 2. Fix for files that DO NOT execute "import torch" but use the word anyway
    setattr(builtins, "torch", mock_instance)
    setattr(builtins, "F", mock_instance)
