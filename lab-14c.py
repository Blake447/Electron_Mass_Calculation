import math

# Global calibration constants
theta_c = 15.250
wavelength_c = 435.833e-9

def main():	
	accepted_mass = 9.10938356e-31 	# accepted mass of the electron

	# initialize the angles used
	violet_blake_left = 164.83
	violet_blake_right = 195.22
	
	violet_hannah_left = 164.85
	violet_hannah_right = 195.20
	
	blue_blake_left = 162.97
	blue_blake_right = 197.10
	
	blue_hannah_left = 162.97
	blue_hannah_right = 197.10
	
	# Calculate the average (unsigned) angle from the center based on left and right inputs
	violet_blake = bisect_angles(violet_blake_left, violet_blake_right)
	violet_hannah = bisect_angles(violet_hannah_left, violet_hannah_right)
	
	blue_blake = bisect_angles(blue_blake_left, blue_blake_right)
	blue_hannah = bisect_angles(blue_hannah_left, blue_hannah_right)
	
	# Average all lab members recorded average angles
	blue_angle = average(blue_blake, blue_hannah)
	violet_angle = average(violet_blake, violet_hannah)
	
	# Calculate the mass of the electron based off the blue and violet angles of diffraction
	electron_mass = electron_mass_from_angles(blue_angle, violet_angle)
	print("mass of electron: " + str(electron_mass) + " kg")
	
	print
	print("Error analysis")
	print
	
	delta = 0.0011
	
	print("Change in blue angle: " +str(delta) + " = largest disagreement among group members (radians)")
	electron_mass_1 = electron_mass_from_angles(blue_angle+delta, violet_angle)
	d_mass_1 = (electron_mass_1 - electron_mass)
	print("Change in mass: " + str(d_mass_1))
	print
	
	print("change in violet angle: " + str(delta))
	electron_mass_2 = electron_mass_from_angles(blue_angle, violet_angle+delta)
	d_mass_2 = (electron_mass_2 - electron_mass)
	print("Change in mass: " + str(d_mass_2))
	print
	
	
	error = pythagorean_sum(d_mass_1, d_mass_2)
	print("Error in mass (Pythogorean sum of change in masses for each wavelength): " + str(error))
	print
	
	print("Results")
	print
	print("Mass of electron: " + str(electron_mass))
	print("Error bars: 	 " + str(error))
	print("accepted mass:	 " + str(accepted_mass))
	print("exper. - theor:   " + str(electron_mass - accepted_mass))
	print("error bars away:  " + str((electron_mass - accepted_mass)/(error)))
	
	print

# Calculate and display information about the electron based on given angles of diffraction from emitted photons
def electron_mass_from_angles(blue_theta, violet_theta):
	print("angle for violet: " + str(violet_theta) + " degrees")
	print("angle for blue:   " + str(blue_theta) + " degrees")

	violet_wavelength = wavelength_from_theta(violet_theta)
	blue_wavelength = wavelength_from_theta(blue_theta)

	print("wavelength for violet: " + str(violet_wavelength) + " m")
	print("wavelength for blue:   " + str(blue_wavelength) + " m")

	electron_mass_violet = electron_mass_from_wavelength(violet_wavelength, 5.0)
	electron_mass_blue = electron_mass_from_wavelength(blue_wavelength, 4.0)
	electron_mass = average(electron_mass_violet, electron_mass_blue)

	electron_mass = apply_correction_factor(electron_mass)	
	return electron_mass

# calculate the mass of the electron from the wavelength of the photon emitted
def electron_mass_from_wavelength(wavelength, level):
	e = 1.6022e-19			# quantum of charge e
	k = 8.9876e9			# coulumbs constant k
	c = 2.9979e8 			# speed of massless particle c
	pi = 3.14159 			# pi
	h = 6.6261e-34/(2.0*pi) 	# reduced plancks constant


	return (4.0*pi*h*h*h*c) / (wavelength*k*k*e*e*e*e*(.25-1.0/(level*level)))

# calculate the wavelength of a diffracted photon given its angle of diffraction
def wavelength_from_theta(theta):
	pi = 3.14159			# circle constant pi
	return (math.sin(theta*2.0*pi/360.0) / math.sin(theta_c*2.0*pi/360.0))*wavelength_c


# apply correction factor
def apply_correction_factor(electron_mass):
	proton_mass = 1.6726e-27
	return electron_mass /(1 - electron_mass / proton_mass)
	
# Find half the angle between two given angles
def bisect_angles(theta_left, theta_right):
	return (theta_right - theta_left)*0.5

# Average two values together
def average(a, b):
	return (a+b)*0.5

# Pythagorean sum
def pythagorean_sum(a, b):
	return math.sqrt(a*a+b*b)

# Run the main program
main()
