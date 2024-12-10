#Kinetic Parameters
kf_C3h2o = 8.3e-7 * 60
kf_C3h2oB = 21.3e4  * 1e-6* 60
kb_C3h2oB = 15.5e-2 * 60
kf_C3h2oH = 5.2e6  * 1e-6* 60
kb_C3h2oH = 32.5 * 60
kb_C3h2oBb = 9.0e-3 * 60
kf_C3bB = 21.3e4  * 1e-6* 60
kb_C3bB = 15.5e-2 * 60
kf_C3Bb = 21.3e4  * 1e-6* 60
kb_C3Bb = 15.5e-2 * 60
kb_C3bBb = 7.7e-3 * 60
kb_C3bBbP = 7.7e-4 * 60
kf_C3bP = 3.0e6  * 1e-6* 60

kb_C3bP = 5.0e-4 * 60
kf_C3bP =3.0e6  * 1e-6* 60
kb_C3bP =5.0e-4 * 60
kf_C3Bsu = 4.2e8  * 1e-6* 60
kf_fC3b = 4.2e8  * 1e-6* 60
kf_hC3b = 4.2e8  * 1e-6* 60
kf_pC3b = 4.2e8  * 1e-6* 60
kf_p⁺re = 1.0e-3 * 60
kf_p⁺su = 3.0e6  * 1e-6* 60
kb_p⁺su =5.0e-4 * 60
kf_iC3bP =3.0e6  * 1e-6* 60
kb_iC3bP =3.8e-4 * 60
kf_C3bH =5.2e6 * 1e-6* 60

kb_C3bH =32.5 * 60
kf_C3bH =5.2e6 * 1e-6* 60
kb_C3bH_ho =3.25 * 60
kf_C3bCR1 =1.2e4 * 1e-6* 60
kb_C3bCR1 =1.0e-2 * 60
kf_C3h2oCR1 =1.2e4 * 1e-6* 60
kb_C3h2oCR1 =1.0e-2 * 60
kf_C3bBbDAF =2.0e3 * 1e-6* 60
kb_C3bBbDAF =7.7e-2* 60
kb_C3bBbCR1 =7.7e-2 * 60
kb_C3bBbH = 7.7e-2 * 60

kf_iC3bCR1 = 2.0e+3* 1e-6* 60
kb_iC3bCR1 =1.0e-2 * 60
kf_C3bBbC3b = 3.5e6 * 1e-6* 60

kb_C3bBbC3b =3.8e-3 * 60
kf_C3bBbC3bC5 = 5.0e6 * 1e-6* 60  
kb_C3bBbC3bC5 =1.0e-2 * 60
kb_C5b =3.8e-2 * 60
kf_C3bBbC3bC5bC6 = 6.0e4 * 1e-6* 60
kb_C3bBbC3bC5bC6 =9e-8 * 60
kf_C5b7 =7.3e5 * 1e-6* 1e-6* 60
kb_C5b7 =1.5e-6 * 60
kf_C5b7su =4.2e8  * 1e-6* 60

kf_mi = 69.3 * 60
kf_C5b8 = 1.1e6 * 1e-6* 60
kb_C5b8 =9.8e-7 * 60
kf_C5b9 = 2.8e6 * 1e-6* 60
kb_C5b9 =2.8e-6 * 60
kf_CnC5b7 =4.1e5  * 1e-6* 60
kb_CnC5b7 = 4.0e-3 * 60
kf_CnC5b8 =4.1e5 * 1e-6* 60
kb_CnC5b8 =4.0e-3 * 60
kf_VnC5b7 =2.4e5  * 60
kb_VnC5b7 =2.0e-3 * 60
kf_CD59C5b9 =1.0e6  * 1e-6* 60
kb_CD59C5b9 =2.0e-4 * 60

kf_CnC5b9 = 4.1e5* 1e-6* 60
kb_CnC5b9 = 4.0e-3* 60
kf_VnC5b8 = 2.4e5* 1e-6* 60
kb_VnC5b8 = 4.0e-3* 60
kf_VnC5b9 = 2.4e5* 1e-6* 60
kb_VnC5b9 = 4.0e-3* 60

kcat_C3h2oBb = 1.8 * 60
km_C3h2oBb = 5.9e-6 *1e6
kcat_C3bBb = 1.8 * 60
km_C3bBb = 5.9e-6 *1e6
kcat_C3bBbP = 3.1 * 60
km_C3bBbP = 1.8e-6 *1e6
kcat_C3bB = 2.1 * 60
km_C3bB = 0.1e-6 *1e6
kcat_C3h2oB = 2.1 * 60
km_C3h2oB = 0.1e-6 *1e6
kcat_C3bH = 1.3 * 60
km_C3bH = 2.5e-7 *1e6
kcat_C3bBbC3b = 4.8 * 60
km_C3bBbC3b = 1.8e-6 *1e6

# kf_C3bBbC3b = 3.5e6 * 1e-6* 60 * 3
# kcat_C3bBbP = 3.1 * 60 * 10

Na = 6.022e23 #Avagadro Number
conc_np = 2e12/Na/1e-3 #Concentration of Nanoparticles (M)

#Diffusion Coefficient Calculation
rc3b = 3.7e-9
Kb = 1.38e-23
T = 310
μ = 4e-3
D = Kb * T / (6 * π * μ * rc3b)

#Local Concentration calculation
t_half = 2e-4
r_hem = √(6*D*t_half) #Radius of 
V_hem = 0.5*4*π*r_hem^3*1000/3 #L
A_hem = π*r_hem^2
A_c3b = 1e-16
c3b_local = A_hem / A_c3b
c_c3b = c3b_local / Na / V_hem
