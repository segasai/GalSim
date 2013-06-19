"""@file test_interpolants_parametric.py  Tests of interpolants using parametric galaxy models.

A companion script to `test_interpolants.py`, but instead of using `RealGalaxy` objects we instead
use Sersic models drawn into `InterpolatedImage` instances to try and get to the nitty-gritty of
the issues with interpolators.

The parameters of the Sersic images come from a COSMOS best-fitting Sersic model catalog.
"""

import cPickle
import numpy as np
import galsim
import test_interpolants


SERSIC_IMAGE_SIZE = 512 # For initial image of the Sersic at Hubble resolution, make nice and large
TEST_IMAGE_SIZE = SERSIC_IMAGE_SIZE  # For speed could make this smaller
# Dictionary for parsing the test_interpolants.interpolant_list into galsim Interpolants 
INTERPOLANT_DICT = {
    "nearest" : galsim.Nearest(),
    "sinc" : galsim.Sinc(),
    "linear" : galsim.Linear(),
    "cubic" : galsim.Cubic(),
    "quintic" : galsim.Quintic(), 
    "lanczos3" : galsim.Lanczos(3),
    "lanczos4" : galsim.Lanczos(4),
    "lanczos5" : galsim.Lanczos(5),
    "lanczos7" : galsim.Lanczos(7)}

# Output filenames
SPACE_FILENAME = 'interpolant_test_parametric_output_space.dat'
GROUND_FILENAME = 'interpolant_test_parametric_output_ground.dat'


class InterpolationDataNoConfig:
    """Quick container class for passing around data from these tests, but not using config.
    """ 
    def __init__(self, g1obs=None, g2obs=None, sigmaobs=None, err_g1obs=None, err_g2obs=None, 
                 err_sigmaobs=None, dx_input=0.03, dx_test=None, shear=None, magnification=None,
                 angle=None, shift=None, x_interpolant=None, k_interpolant=None, pad_factor=None,
                 image_type=None, psf_label='space'):
        self.g1obs = g1obs
        self.g2obs = g2obs
        self.sigmaobs = sigmaobs
        self.err_g1obs = err_g1obs
        self.err_g2obs = err_g2obs
        self.err_sigmaobs = err_sigmaobs
        self.dx_input = dx_input
        self.dx_test = dx_test
        self.test_image_size = TEST_IMAGE_SIZE
        # Parse and store default/non-default shear, mag, rotation and shift
        if shear is None:
            self.shear = [0., 0.]
        else:
            self.shear = shear
        if magnification is None:
            self.magnification = 1.
        else:
            self.magnification = magnification
        if angle is None:
            self.angle = 0.
        else:
            self.angle = angle.rad * 180. / np.pi
        if shift is None:
            self.shiftx = 0.
            self.shifty = 0 .
        else:
            self.shiftx = shift.x
            self.shifty = shift.y
        # Store the interpolants
        if x_interpolant is None:
            self.x_interpolant = 'default'
        else:  
            self.x_interpolant = x_interpolant
        if k_interpolant is None:
            self.k_interpolant = 'default'
        else:  
            self.k_interpolant = k_interpolant
        # Store padding & image type
        if pad_factor is None:
            self.padding = 0
        else:
            self.padding = pad_factor
        self.image_type = "Parametric (space OR ground)"
        self.psf_label = psf_label


def calculate_interpolated_image_g1g2sigma(images, psf=None, dx_input=None, dx_test=None, 
                                           shear=None, magnification=None, angle=None, shift=None, 
                                           x_interpolant=None, k_interpolant=None, padding=None):
    """Takes a list of drawn images of Sersic galaxies, reads them into an InterpolatedImage object
    using the supplied parameters, and calculates the g1, g2 and sigma of the output.
    """
    # Some input parsing
    if padding is None:
        pad_factor = 0.
    else:
        pad_factor = padding
    # Loop over the images and generate an InterpolatedImage from the pixel values
    g1obs_list = []
    g2obs_list = []
    sigmaobs_list = []
    for sersic_image in sersic_images:

        # Build the raw InterpolatedImage
        test_gal = galsim.InterpolatedImage(
            sersic_image,
            dx=dx_input,
            x_interpolant=INTERPOLANT_DICT[x_interpolant],
            k_interpolant=INTERPOLANT_DICT[k_interpolant],
            pad_factor=pad_factor)
        # Apply shears, magnification, rotation and shifts if requested
        if shear is not None:
            test_gal.applyShear(g1=shear[0], g2=shear[1])
        if magnification is not None:
            test_gal.applyMagnification(magnification)
        if angle is not None:
            if not isinstance(angle, galsim.Angle):
                raise ValueError("Input kwarg angle must be a galsim.Angle instance.")
            test_gal.applyRotation(angle)
        if shift is not None:
            if not isinstance(shift, galsim.PositionD):
                raise ValueError("Input kwarg shift must be a galsim.PositionD instance.")
            test_gal.applyShift(shift) 
        # Apply a PSF if requested
        if psf is not None:
            test_final = galsim.Convolve([test_gal, psf])
        else:
            test_final = test_gal
        # Draw into the test image and calculate adaptive moments
        test_image = galsim.ImageD(TEST_IMAGE_SIZE, TEST_IMAGE_SIZE)
        test_image.setScale(dx_test)
        test_final.draw(test_image, dx=dx_test)
        trial_result = test_interpolants.CatchAdaptiveMomErrors(test_final)
        if isinstance(trial_result, float):
            g1obs_list.append(-10)
            g2obs_list.append(-10)
            sigmaobs_list.append(-10)
        elif isinstance(trial_result, galsim.hsm.ShapeData):
            g1obs_list.append(trial_result.observed_shape.g1)
            g2obs_list.append(trial_result.observed_shape.g2)
            sigmaobs_list.append(trial_result.moments_sigma)
        else:
            raise TypeError("Unexpected output from test_interpolants.CatchAdaptiveMomErrors().")

    ret = InterpolationDataNoConfig(
        g1obs=g1obs_list, g2obs=g2obs_list, sigmaobs=sigmaobs_list, dx_input=dx_input,
        dx_test=dx_test, shear=shear, magnification=magnification, angle=angle, shift=shift,  
        x_interpolant=x_interpolant, k_interpolant=k_interpolant, padding=padding)
    return ret

def draw_sersic_images(narr, hlrarr, gobsarray, random_seed=None, nmin=0.3, nmax=4.2,
                       image_size=512, pixel_scale=0.03):
    """Given input NumPy arrays of Sersic n, half light radius, and |g|, draw a list of Sersic
    images with n values within range, at random orientations.
    """
    # Initialize the random number generator
    if random_seed is None:
        u = galsim.UniformDeviate()
    else:
        u = galsim.UniformDeviate(random_seed)

    # Loop over all the input params and make an sersic galaxy image from each
    sersic_images = []
    print "Drawing Sersic images"
    for n, hlr, gobs in zip(narr, hlrarr, gobsarr):

        # First check we are only using Sersic indices in range
        if n <= nmin or n >= nmax: continue # This goes to the next iteration

        # Otherwise set up the image to draw our COSMOS sersic profiles into
        sersic_image = galsim.ImageD(image_size, image_size)
        sersic_image.setScale(pixel_scale)
        # Build the galaxy
        galaxy = galsim.Sersic(n=n, half_light_radius=hlr)
        # Apply the ellipticity of the correct magnitude with a random rotation
        theta_rot = 2. * np.pi * u() # Random orientation
        galaxy.applyShear(g1=gobs*np.cos(2.*theta_rot), g2=gobs*np.sin(2.*theta_rot))
        galaxy.draw(sersic_image, dx=test_interpolants.space_pixel_scale)
        sersic_images.append(sersic_image)

    # Return this list of drawn images
    return sersic_images


if __name__ == "__main__":

    # Import the Sersic galaxy sample module
    try:
        import galaxy_sample
    except ImportError:
        import sys
        sys.path.append('../external/test_sersic_highn')
        import galaxy_sample

    # Get the COSMOS galaxy sample parameters
    ns_cosmos, hlrs_cosmos, gobss_cosmos = galaxy_sample.get()
    # Only use the first test_interpolants.nitems galaxies in these lists, starting at
    # test_interpolants.first_index
    istart = test_interpolants.first_index
    iend = istart + test_interpolants.nitems
    ns_cosmos = ns_cosmos[istart: iend]
    hlrs_cosmos = hlrs_cosmos[istart: iend]
    gobss_cosmos = gobss_cosmos[istart: iend]

    # Draw a whole load of images of Sersic profiles at random orientations using these params
    sersic_images = draw_sersic_images(
        ns_cosmos, hlrs_cosmos, gobss_cosmos, random_seed=test_interpolants.rseed, nmin=0.3,
        nmax=4.2, image_size=SERSIC_IMAGE_SIZE, pixel_scale=test_interpolants.space_pixel_scale)

    # Let's just do space and ground-based PSFs, and define a dict storing these to iterate over,
    # along with the appropriate test pixel scale
    psf_dict = {
        "space" : (
            galsim.Airy(lam_over_diam=test_interpolants.space_lam_over_diam),
            test_interpolants.space_pixel_scale),
        "ground" : (
            galsim.Convolve(
                galsim.Kolmogorov(fwhm=test_interpolants.ground_fwhm),
                galsim.Airy(lam_over_diam=test_interpolants.ground_lam_over_diam)),
            test_interpolants.ground_pixel_scale)
    }

    # Then we start the grand loop producing output in a similar fashion to test_interpolants.py
    ground_file = open(GROUND_FILENAME, 'wb')
    space_file = open(SPACE_FILENAME, 'wb')
    for psf_label in ("space", "ground"):
 
        # Get the correct PSF and test image pixel scale
        psf = psf_dict[psf_label][0]
        dx_test = psf_dict[psf_label][1]

        for padding in test_interpolants.padding_list:

            for interpolant in test_interpolants.use_interpolants:

                print 'Calculating unrotated, undistorted observed shears and sizes'
                intdata0x = calculate_interpolated_image_g1g2sigma(
                    sersic_images, psf=psf, dx_input=test_interpolants.space_pixel_scale,
                    dx_test=dx_test, shear=None, magnification=None, angle=None, shift=None, 
                    x_interpolant=interpolant, padding=None)
                intdata0k = calculate_interpolated_image_g1g2sigma(
                    sersic_images, psf=psf, dx_input=test_interpolants.space_pixel_scale,
                    dx_test=dx_test, shear=None, magnification=None, angle=None, shift=None, 
                    x_interpolant=interpolant, padding=None)
                
