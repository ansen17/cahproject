Determining Contact Angle Hysteresis on Sinusoidal Physically Patterned Surfaces
==============

Summary of the Project
--------------

I completed a computational modeling project for my undergraduate honors thesis project. I focused on determining a thermodynamic model for droplets
on sinusoidal physically patterned surfaces. I used Python scripts to help obtain relevant data about the droplets in conjunction with a C-based 
program called *Surface Evolver*. I developed an droplet evolution algorithm and various phase diagrams to show how energy related to the droplet 
behavior. Below is a broad description of what I did. 

I have hosted most of the core files relevant in obtaining the final result between energy and contact angle hysteresis (a physical parameter). Other 
files not included consisted largely of automation and data retrieval Python scripts.

Please let me know if you have any questions and thank you for checking my page.



Background of Wetting

Wetting describes a fluid's ability to remain contact with a surface and contact angle and its hysteresis are the most common quantitative measures 
to illustrate this phenomenon. Contact angle hysteresis is defined as the receding contact angle subtracted from the advancing contact angle. Contact 
angle hysteresis’ significance is that for low contact angles allow for droplets of liquids to slide off surfaces very easily and similarly for large 
contact angle hysteresis, it is harder for liquid droplets to slide off surfaces amounting to a variety of industrial and research applications. Various 
characteristics of a sinusoidal patterned surface like the wavelength, amplitude, grooves and energetic properties like the strength of gravity and type 
of thermodynamic wetting model affect contact angle hysteresis of the droplet and therefore the sliding behavior of the droplet. This work can be used as 
a basis for other patterned surfaces like triangular, rectangular, and hierarchal rough surfaces in the future. The importance of contact angle hysteresis 
and sliding droplets is that the research development in this field can help produce valuable contributions in designing a variety of advanced and soft 
materials that have applications in a multitude of fields from microfluidics to coating.

*Specifics of Project*

- Edited the evolver.py file that initiated the automation process of sending thousands of files to LionX cluster
- Added additional functions to .fe file which is the file that the program *Surface Evolver* (a C-based program) that depicts the droplet and the 
  energy functions acting on it
- After trials and data procurement, established a set of commands in Python script that obtained several measurements at each tilt angle of the plane from 0 
  to 90 degrees
- Data measurements included amplitude of the surface, wavelength of the surface, bond number, tilt angle of the plane, advancing and receding contact angles, 
  contact angle hysteresis and the unitless scaled energy value that *Surface Evolver* would output 
- Create phase diagrams between different sets of data sets

   