class Config:
    """
    Config class to hold parameters for the biexponential transformation.

    Attributes:
    ----------
    width_multi : int
        The width of the linear region of the transformation. This parameter controls the 
        smoothness of the transition between the linear and logarithmic regions.

    pos : float
        The positive minimum value for the transformation. This value defines the lower 
        bound for positive data in the transformation, essentially controlling where the 
        transition to the linear region starts. Default is 4.9.

    neg : int
        The negative minimum value for the transformation. This sets the minimum threshold 
        for the negative side of the data, which is important for data that might include 
        negative values or for symmetric scaling. Default is 0.
    """
    width_multi:int = 0.00001
    pos:float = 4.418540 
    neg:float = 0.5 # range 0-1