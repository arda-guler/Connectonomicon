=== === === CONNECTONOMICON README === === ===
version 20241219.00

!!!!!
If you don't want to read this, please just 
take a look at "Important Notes" at line 206.
!!!!!

=== === Description
Connectonomicon is a tool created to assist 
with the manual identification of minor planet 
observation tracklets.

=== === Requirements
Connectonomicon is a Python 3 script, but it
relies on the Windows console version of 
find_orb. However, although it is not
recommended, you can use the native orbit
determination and run it on any OS that
can run Python 3.

Connectonomicon makes use of the following 
3rd party Python packages:

numpy (tested with 1.24.2)
spiceypy (tested with 6.0.0)
matplotlib (tested with 3.5.3)
mplcursors (tested with 0.6)

(+ some native Python modules in Python 3.10)

=== === Setting Up
This is likely the most boring part.

1) find_orb (Optional, highly recommended)
It is recommended that you download the console
version of find_orb (fo64.exe) from the
following URL:

https://www.projectpluto.com/fo_usage.htm
(The link that reads 'click here for a 
pre-built fo64 for 64-bit Windows'.)

Place find_orb in the same directory as
Connectonomicon.

2) Obtain SPICE kernels
Obtain the following from NAIF:
de440.bsp
naif0012.tls
earth_000101_250316_241218.bpc
pck00011.tpc

https://naif.jpl.nasa.gov/naif/data.html

Place them in a 'data' folder:
data/de440.bsp
data/naif0012.tls
data/earth_000101_250316_241218.bpc
data/pck00011.tpc

3) Obtain observatory data from MPC
Save the following page

https://minorplanetcenter.net/iau/lists/ObsCodes.html

as data/obscodes.txt.

4) Install required Python packages
Using pip is the most common and easiest way.
Go to the command line and use:

pip install numpy
pip install spiceypy
pip install matplotlib
pip install mplcursors

...or, if you need to install specific versions:

pip install numpy==1.24.2
pip install spiceypy==6.0.0
pip install matplotlib==3.5.3
pip install mplcursors==0.6

See this if you are stuck:
https://docs.python.org/3/installing/index.html

Now, you are set!

=== === How to Use
1) Preparing input files
You need two files holding astrometry data in
Obs80 format (MPC's 80-column format). One of
them will hold the known observations of the
object being investigated, while the other will
hold the list of observations you wish to
classify and link. (eg. a part of MPC's ITF file)
You can name the files as you wish.

2) Running the program
Run the script to launch a small GUI screen.
Enter the inputs described below:

Known Observations File: name of the file
    holding the observations belonging to the
    object currently being investigated.
    
Potential Observations File: name of the file
    holding unlinked astrometric measurements
    you wish to classify
    
Propagation Time: number of days you wish to
    simulate the estimated orbit of the
    object for, after the final known observation
    
Max. Measurement Error: max. possible measurement 
    error for the object's position on the image
    when the astrometry was created
    
Pixel Resolution: Angular resolution of the
    electrooptical system through which the
    measurements were made (often a CCD cam).
    
Reference Obscode: The MPC obscode of the
    observatory for which you want to plot the
    chart. It won't make much difference for
    distant objects but may be especially 
    important for NEOs and is practically a must 
    for artsats.
    
Orbit Determination: Orbit determination method.
    find_orb is strongly recommended.
    
Plot Unlikely Observations: It is suggested to
    keep this unchecked if the number of your 
    unlinked observations are too high as it 
    impacts performance.
    
3) Run analysis
Click on the button that reads 'RUN' and wait.
The program will inform you of its progress via
console printouts. It should eventually create
an interactive plot.

4) Examine the plot
Use the controls on the plot screen to zoom and
pan the view. Click on individual observations
to see their 80-column formats, and try fitting
orbits with the tracklets that you think may
belong to your object.

== == How Does It Work?
For a slightly more detailed explanation, see
Connectonomicon.pdf.

The Obs80 format uses fixed locations for various
data such as observation time, RA, DEC eg., which
makes it very easy to parse them. The program
internally uses degrees for angles and 'datetime'
objects for times.

If find_orb is used as the orbit determinator:
The file that holds the known observations is
loaded into find_orb, which silently generates
an elements.txt file. In that file is also the
estimated state vectors at some epoch. 
Connectonomicon reads it for use in its orbit
propagator.

If native built-in orbit determinator is used:
The first and last observation, along with
another observation towards the middle of the
list is taken. The program initially uses the
Laplace method until it obtains an initial
distance estimate for the object (from the
observer). Then it switches to a modified
Herget method where it gradually nudges the
position and velocity vectors to reduce the
RA-DEC errors as much as possible. The
estimated initial state vector is then used
for orbit propagations.

The orbit propagation simulates the orbit of
the object. The force model includes the Sun
and the planet-moon system barycenters. The 
propagator is an 8th order leapfrog-ish method 
(Yoshida-type). Non-gravitational perturbations
are not accounted for, because the error one
will get due to the state vector estimation
will easily make non-gravs irrelevant (unless
you are looking at an alien spaceship).

To aid with short-arc estimations, a
straight-line constant-sky-speed estimation
is also done. Based on your max. pixel error
estimations, an area of the sky will be
determined where the next observation may be
found after your given number of propagation 
days.

After orbit propagations, the unlinked
observations will be classified according to
their dates, their proximity to the estimated
orbit path, and their magnitudes.

== == Important Notes!   
 > find_orb handles observation sets with
   different temporary codes as different
   objects. Try using the same temporary code
   for all known observations.
   
 > Orbit propagation errors will increase as
   propagation time is increased. Chances are,
   if you are interested in tracklet IDing,
   you don't need me to tell you this.
   
 > If the object's magnitude changes rapidly,
   actual observations of the object may get
   filtered out. Mark 'Plot Unlikely
   Observations' while working with such 
   objects.
   
 > The native orbit determination algorithm
   is experimental!
   
== == Contact
Contact the author via email, or open an issue on
Connectonomicon GitHub repository. You can
contact for bug reports, fixes, feature requests,
or to say hi. 

I do astronomy in my free time, so I might take 
a bit to reply back.

a r d a g u l e r 0 9 (gmail)
arda-guler (GitHub)

== == Connectonomicon License
Connectonomicon
Copyright (C) 2024  H. A. Guler

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see
<https://www.gnu.org/licenses/>.

== == Acknowledgements
Connectonomicon makes use of find_orb orbit
determination software by Bill Gray. (Thanks!)

https://www.projectpluto.com/find_orb.htm

Connectonomicon uses SPICE kernels provided by
NAIF as part of its orbit propagation scheme:

https://naif.jpl.nasa.gov/naif/credit.html

You should have received neither find_orb
nor any SPICE kernels with Connectonomicon.
