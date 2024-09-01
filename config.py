class Config:
    """
    Config class to hold parameters for the biexponential transformation.

    Attributes:
    ----------
    max_value : int
        The maximum value for the data to be transformed. This value sets the upper limit 
        for the transformation. Typically, this is set to a large number to encompass the 
        entire range of your data. Default is 10,000,000.

    width : int
        The width of the linear region of the transformation. This parameter controls the 
        smoothness of the transition between the linear and logarithmic regions. A negative 
        value can be used to adjust the width in a way that suits the data characteristics. 
        Default is -100.

    pos : float
        The positive minimum value for the transformation. This value defines the lower 
        bound for positive data in the transformation, essentially controlling where the 
        transition to the linear region starts. Default is 4.9.

    neg : int
        The negative minimum value for the transformation. This sets the minimum threshold 
        for the negative side of the data, which is important for data that might include 
        negative values or for symmetric scaling. Default is 0.
    """
    max_value:int = 10000000
    width:int = -100
    pos:float = 4.9 
    neg:float = 0 # range 0-1