# DRASTIC
**D**une c**R**p **A**node**S** pho**T**o qual**I**ty **C**ontrol 

Julia software developed at LAPP to perform quality assurance/control of DUNE CRP anodes based on image analysis

# Usage
To use the software you need ot make available in you workspace Julia, which is available at lxplus, cc in2p3, etc 
Once you have make Julia available, you can simple install the DRASTIC package by doing inside julia
```
]add https://github.com/lmanzanillas/DRASTIC
```
The previous command will install the DRASTIC package 
In addition to the DRASTIC package we recomend to install the following complementary packages for analyzing the results
```
]add Images, ImageFeatures, FileIO, Plots, ImageComponentAnalysis, Statistics, StatsBase, HDF5
```

# Step by step example
The Quality control of the CRP anodes will check if the pitch of the PCB holes correspond to the required specifications and if no anomalous pattern is observed
To perform this QC photos will be taken of the PCBs/anodes that will be analyzed to extract the required paramters for the analysis
Since the axes of the 2D images are in pixels we will need to find a calibration factor to convert pixel units to mm. To this end a photo with calibration red circles will be used. The diameter of this circles is about 10 mm.

Here photo

Then we need to define the color that will be used as reference for the calibration, i.e. the red color. To this end, we will select a region of the red circles that we will convert to the HSV color space. The photo is a matrix of 3456 x 5184 pixels. You can use ```plot(c_img[1900:2050,4300:4500])``` to make sure are selecting the good region.

Here red color

You can use the same image to define the hole color, which corresponds to the center of the holes (the blueish part). The idea will be take a hole placed in the center of the image. In this case: 
```python
hole_color = HSV{Float32}( mean(c_img[1505:1550,2635:2685]))
```

## Calibration factor
The first parameter that we need to find is the calibration factor. To this end the function
```
get_calibration_factor ```
can be used. This function need an image with the red circles and a color as reference, in this case the red color
```python
calib = get_calibration_factor(c_img,red_color)
```

## Analysis  of a photo
Once the calibration factor has been found we can proceed with the analysis
First we read a photo for the analysis
```python
light_dir = "/home/manzanilla/Pictures/DUNE_PCBs_QA/50mm_test/"
my_photos = light_dir .* filter(x->occursin("002.jpg",x), readdir(light_dir))
testing_img = load(my_photos[1])
```

here photo

Since the light can be different in different regions of the photo, we will divide the photo in small section to have more homogenous results in all the photo. To do that, you can use the following function
```python
div_img, indexv = divide_img_sq(testing_img,430);
```
which will divide the photo in section of 430x430 pixels (except the last section that will be adjusted to keep the total amount of pixels)
This function will return a vector containing the sections and a matrix containing the position of each section to allow the merging of all sections into the original image


