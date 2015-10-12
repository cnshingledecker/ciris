PRO cross_sections

; This section calculates the ionization cross-section
A = 2.98
B = 4.42
C = 1.48
D = 0.750
F = 4.80
tau = 100E3 ; the proton kinetic energy in eV (100 keV)
bohr = 5.2918E-11 ; the Bohr radius in m
ryd = 1.0973731E7 ; the Rydberg constant in m^-1
pi = 3.14159 ; Pi
m_e = 9.10938E-31 ; electron mass in kg
m_p = 1.672622E-27 ; proton mass in kg
T = (m_e/m_p)*tau ; the kinetic energy of an electron equivalent to a proton

; Calculate the low-energy cross-section
sigma_low = 4*pi*(bohr^2)*(C*((T/ryd)^D)+F)
print, "The low energy cross-section is",sigma_low," m^2"

; Calculate the high-energy cross-section
sigma_high = 4*pi*(bohr^2)*(ryd/T)*(A*ALOG(1+(ryd/T))+B)
print, 'The high energy cross-section is',sigma_high,' m^2'

; Calculate the total cross-section, which is the harmonic mean of the high and low
sigma_ion = ((1/sigma_low)+(1/sigma_high))^(-1)
print, 'The total ionization cross-section is', sigma_ion,' m^2'


; This section calculates the excitation cross-sections and plots
; them as a function of # of electrons in the target material
C1 = 4
C2 = 0.25
Z = 16; the number of electrons in the target material (O2)
Ek = 7.083 ; the excitation energy leading to dissociation in oxygen from Cosby 1993 in eV
sigma_0 = 10E-20 ; a constant in m^2
nu = 1 ; from Miller and Green (1973) via Dingfelder 
omega = 0.85 ; same as above
sigma_e = 5E-22 ; 

J = C2*(m_p/m_e)*T*(omega/nu)^(1/(omega+nu))
a = C2*(m_p/m_e)*(T/Z)*((C1)/(sigma_0*nu))^(1/omega)


END 