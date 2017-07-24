# curvefit
A command line tool for nonlinear regression analysis using the
Levenberg-Marquardt algorithm.

## Supported Models
    - Boltzmann sigmoid:

      y = min + ((max - min)/(1 + exp((v50 - x)/slope)))

    - Exponential decay:

      y = a + b * exp(-lambda * x)

    - Gaussian function:

      y = a + b * exp(-mu^2/sigma^2)

    - Hill equation (allosteric kinetics)

      y = delta / (1 + (ic50 / x)^hill)

    - Dose response (ic50):

      y = 1 - (delta / (1 + (ic50 / x)^hill))

    - Michaelis-Menten:

      y = Vmax * x / (Km + x)

    - Modified sine wave

      y = a * sin(pi * (x - b) / c)

The program finds values for the parameters of the model function
(e.g. a, b, c, d) that minimizes the difference between input and 
fitted values of y for all values of x.

## Requirements
[GNU Scientific Library](http://www.gnu.org/software/gsl/)


## Usage
The command to run the program is

`./curvefit -f datafile -m <model> [-a <float> -b <float> -c <float> -d <float>]`

## Required Arguments
-f      Path to input data. Must be a tab delimited text file of (x, y) value 
        pairs. No other data should be included in the file.

-m      Desired model. Valid options are: 

        boltzmann   -   Boltzmann sigmoid
        expdecay    -   Exponential decay
        gaussian    -   Gaussian
        hill        -   Hill equation
        ic50        -   Dose Response
        mm          -   Michaelis-Menten
        modsin      -   Sine wave

## Optional Arguments
-abcd   Specify initial values for the parameters to fit. Default values are
        1.0.

        flag    boltzmann expdecay  gaussian  hill     ic50     mm      modsin
        ----------------------------------------------------------------------
        -a      min       a         a         delta    delta    Vmax    a
        -b      max       b         b         ic50     ic50     Km      b
        -c      v50       lambda    mu        hill     hill     N/A     c
        -d      slope     N/A       sigma     N/A      N/A      N/A     N/A
