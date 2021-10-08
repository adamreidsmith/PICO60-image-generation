# Generating Artificial Event Images from the PICO-60 Dark Matter Experiment

This code analyzes images and generates artificial optical data from the [PICO-60 dark matter search experiment](https://en.wikipedia.org/wiki/PICO) using a generative adversarial network (GAN).  PICO-60 is a [bubble chamber](https://en.wikipedia.org/wiki/Bubble_chamber) experiment searching for direct evidence of WIMP dark matter particles which was operational in 2016 and 2017.  The detector utilized 4 high speed cameras to image a 60 litre volume of superheated liquid freon with the aim of observing bubbles produced by nuclear recoils from WIMP-nucleon scattering.  Due to its low-background design, PICO-60 lacks the large quantities of data necessary for training of machine learning algorithms.  This project provides a method for generating artificial optical data mimicking images obtained from the PICO-60 detector’s four cameras.

A comprehensive analysis of the project, including the theory, methods, and analysis of results, is provided in [`Artificial_Image_Generation.pdf`](/Artificial_Image_Generation.pdf).

## Usage

A comprehensive overview of the usage of the code in this project is provided in [`image_generation_technote.pdf`](/image_generation_technote.pdf).

## Technologies used

  * Python 3
  * MATLAB
  * PyTorch
  * OpenCV
  * MATLAB Engine API for Python
  * NumPy
  * SciPy
  * Matplotlib

## Credits

This project was produced by Adam Smith for the Ma Ph 499 Undergraduate Research Project course at the University of Alberta under supervision of Dr. Carsten Krauss.  All python code and analysis (with the exception of [`setup.py`](/setup.py)) was written and performed by Adam Smith.  The MATLAB code provided in [`raytracerfiles/`](/raytracerfiles) was originally written by Eric Dahl and later modified by Gavin Crowder and Clarke Hardy.  The files [`FitPICO60Geometry.m`](/raytracerfiles/FitPICO60Geometry.m) and [`GetRaysAndPixels.m`](/raytracerfiles/GetRaysAndPixels.m) were adapted by Adam Smith from [`runthenewgeometryfitter_pico60_2016.m`](/raytracerfiles/runthenewgeometryfitter_pico60_2016.m) and [`EricsQuickLookupTableMakerEtcWithTorus.m`](/raytracerfiles/EricsQuickLookupTableMakerEtcWithTorus.m) respectively for use with this project.

## License

[MIT](LICENSE)
