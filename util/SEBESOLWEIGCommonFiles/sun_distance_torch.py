__author__ = "xlinfr"
try:
    import torch
except:
    pass


def sun_distance(jday):
    """

    #% Calculatesrelative earth sun distance
    #% with day of year as input.
    #% Partridge and Platt, 1975
    """
    b = 2.0 * torch.pi * jday / 365.0
    D = torch.sqrt(
        (
            1.00011
            + 0.034221 * torch.cos(b)
            + 0.001280 * torch.sin(b)
            + 0.000719 * torch.cos(2.0 * b)
            + 0.000077 * torch.sin(2.0 * b)
        )
    )
    return D
