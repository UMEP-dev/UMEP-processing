def Kvikt_veg(input_svf, svfveg, vikttot):

    # Least
    viktwall = (
        vikttot
        - (
            63.227 * input_svf**6
            - 161.51 * input_svf**5
            + 156.91 * input_svf**4
            - 70.424 * input_svf**3
            + 16.773 * input_svf**2
            - 0.4863 * input_svf
        )
    ) / vikttot

    svfvegbu = svfveg + input_svf - 1  # Vegetation plus buildings
    viktveg = (
        vikttot
        - (
            63.227 * svfvegbu**6
            - 161.51 * svfvegbu**5
            + 156.91 * svfvegbu**4
            - 70.424 * svfvegbu**3
            + 16.773 * svfvegbu**2
            - 0.4863 * svfvegbu
        )
    ) / vikttot
    viktveg = viktveg - viktwall

    return viktveg, viktwall
