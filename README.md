# SMALL-LABS

The SMALL-LABS (Single-Molecule Accurate Localization by Local Background Subtraction) algorithm, accurately locates and measures the intensity of single molecules, regardless of the shape or brightness of the background.

The program can also fit single molecules without doing background subtraction. See the User Guide for more details.

Written by Benjamin P Isaacoff at the University of Michigan.

## Installation

Download the entire folder and unzip if you downloaded the .zip folder. Change the working directory in Matlab to this folder and call the functions in the Matlab command window as described in the User Guide.

## Usage

See the Quick Start Guide for a quick introduction to using *SMALLLABS_main*. 

See the User Guide for the details. Briefly, the function *SMALLLABS_main* is a wrapper for the other code to perform all of the steps in the correct order. Simply run *SMALLLABS_main* by specifying the directory containing your movies, specify the three required parameters, and any optional parameters, then run it and click to choose the movies you want to fit. Or run the various programs independently.

## Contributing

Please inform us (isaacoff@umich.edu) before making any changes, then follow the directions below: 
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## Credits

All individual programs should have all their individual attributions still in place including authors. 

This code was developed with support from the National Science Foundation (NSF grant CHE-1252322).

The development of this code is greatly indebted to the work of David J Rowland (often referred to as DJR in the code). In addition to containing some functions written by him, I’ve borrowed a lot of code snippets from his programs.

The code not written by BPI:

*TiffStack* by DR Muir and BM Kampa  
*saveastiff* by YoonOh Tak  
*bpass* by John C. Crocker and David G. Grier  
*MLEwG* by KI Mortensen, LS Churchman, JA Spudich, H Flyvbjerg  
*gaussFit* by David J Rowland  
*Track_3D2* by David J Rowland   
*hungarian* by Yi Cao  
*gpufit* by Adrian Przybylski, Björn Thiel, Jan Keller-Findeisen, Bernd Stock, and Mark Bates

## License

                      GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

  See LICENSE.txt
