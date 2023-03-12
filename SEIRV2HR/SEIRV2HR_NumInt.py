import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class SEIRSV2_Model():

    def __init__(self):

        #General Coefficients
        self.PopulationHR           = 0
        self.PopulationLR           = 1
        self.PopulationTotal        = self.PopulationLR + self.PopulationHR  
        self.InitialExposuresHR     = 0
        self.InitialExposuresLR     = 0.00005
        self.InitialExposures       = self.InitialExposuresHR + self.InitialExposuresLR

        self.VaccineDay             = 0
        self.BirthRateLR            = 2.78363e-5
        self.BirthRateHR            = 0 #2.78363e-5
        self.DeathRate              = 0 #2.48321e-5
        
        #Pandemic Coefficients
        self.InfectionRate          = 0.5
        self.RecoveryRate           = 0.05
        self.ImmunityLossRate       = 0.0035
        self.ContactRateS           = 1/5
        self.ContactRateSV          = 0 #self.ContactRateS * 0.19
        self.ContactRateDV          = 0 #self.ContactRateS * 0.09
        self.ContactRateSHR         = self.ContactRateS  * 2
        self.ContactRateSVHR        = self.ContactRateSHR * 0.19
        self.ContactRateDVHR        = self.ContactRateSHR * 0.09
        self.VaccinationRateS       = 0.01
        self.VaccinationRateSV      = 0.01
        self.VaccinationRateSHR     = 0
        self.VaccinationRateSVHR    = 0
        self.DeathRateCOVID         = 0 #8.95027e-4
        self.DeathRateCOVIDHR       = 0 #self.DeathRateCOVID * 2
        self.DeathRateCOVIDSV       = 0 #self.DeathRateCOVID * 0.19
        self.DeathRateCOVIDDV       = 0 #self.DeathRateCOVID * 0.09
        self.DeathRateCOVIDSVHR     = 0 #self.DeathRateCOVIDHR * 0.19
        self.DeathRateCOVIDDVHR     = 0 #self.DeathRateCOVIDHR * 0.09

    def NumInt(self):

        S = [self.PopulationLR - self.InitialExposuresLR]; E = [self.InitialExposuresLR]; I = [0]; R = [0]
        SV  = [0]; SVE = [0]; SVI = [0]
        DV  = [0]; DVE = [0]; DVI = [0]
        SHR = [self.PopulationHR - self.InitialExposuresHR]; EHR = [self.InitialExposuresHR]; IHR = [0]; RHR = [0]
        SVHR = [0]; SVEHR = [0]; SVIHR = [0]
        DVHR = [0]; DVEHR = [0]; DVIHR = [0]

        for n in np.arange(0, 1000):

            C = I[n] + E[n] #+ SVE[n] + SVI[n] + DVE[n] + DVI[n]
            
            dS = self.BirthRateLR*self.PopulationTotal - self.ContactRateS*(C/self.PopulationTotal)*S[n] + self.ImmunityLossRate*R[n] - self.DeathRate*S[n] - self.VaccinationRateS*S[n]
            dE = self.ContactRateS*(C/self.PopulationTotal)*S[n] - self.InfectionRate*E[n] - self.DeathRate*E[n] 
            dI = self.InfectionRate*E[n] - self.RecoveryRate*I[n] - self.DeathRateCOVID*I[n] 
            dR = self.RecoveryRate*I[n] - self.ImmunityLossRate*R[n] - self.DeathRate*R[n] - self.VaccinationRateS*R[n]

            dSV = self.VaccinationRateS*S[n] + self.VaccinationRateS*R[n] - self.VaccinationRateSV*SV[n]
            #dSVE = self.ContactRateSV*(C/self.PopulationTotal)*SV[n] - self.InfectionRate*SVE[n] - self.DeathRate*SVE[n]
            #dSVI = self.InfectionRate*SVE[n] - self.RecoveryRate*SVI[n] - self.DeathRateCOVIDSV*SVI[n]

            dDV = self.VaccinationRateSV*SV[n] #- self.ContactRateDV*(C/self.PopulationTotal)*DV[n] + self.RecoveryRate*DVI[n] - self.DeathRate*DV[n] 
            #dDVE = self.ContactRateDV*(C/self.PopulationTotal)*DV[n] - self.InfectionRate*DVE[n] - self.DeathRate*DVE[n]
            #dDVI = self.InfectionRate*DVE[n] - self.RecoveryRate*DVI[n] - self.DeathRateCOVIDDV*DVI[n]


            #C = I[n] + E[n] + SVI[n] + SVE[n] + DVI[n] + DVE[n] + EHR[n] + IHR[n] + SVIHR[n] + SVEHR[n] + DVEHR[n] + DVIHR[n] 

            #dS      = self.BirthRateLR*self.PopulationTotal - self.ContactRateS*(C/self.PopulationTotal)*S[n] + self.ImmunityLossRate*R[n] - self.VaccinationRateS*S[n]  - self.DeathRate*S[n] 
            #dE      = self.ContactRateS*(C/self.PopulationTotal)*S[n] - self.InfectionRate*E[n] - self.DeathRate*E[n] 
            #dI      = self.InfectionRate*E[n] - self.RecoveryRate*I[n] - self.DeathRateCOVID*I[n]
            #dR      = self.RecoveryRate*I[n] - self.VaccinationRateS*R[n] - self.ImmunityLossRate*R[n] - self.DeathRate*R[n]

            #dSV     = self.VaccinationRateS*(S[n] + R[n]) - self.ContactRateSV*(C/self.PopulationTotal)*SV[n] - self.DeathRate*SV[n] - self.VaccinationRateSV*SV[n] + self.RecoveryRate*SVI[n]
            #dSVE    = self.ContactRateSV*(C/self.PopulationTotal)*SV[n] - self.InfectionRate*SVE[n] - self.DeathRate*SVE[n]
            #dSVI    = self.InfectionRate*SVE[n] - self.RecoveryRate*SVI[n] - self.DeathRateCOVIDSV*SVI[n]

            #dDV     = self.VaccinationRateSV*SV[n] - self.ContactRateDV*(C/self.PopulationTotal)*DV[n] + self.RecoveryRate*DVI[n] - self.DeathRate*DV[n] 
            #dDVE    = self.ContactRateDV*(C/self.PopulationTotal)*DV[n] - self.InfectionRate*DVE[n] - self.DeathRate*DVE[n]
            #dDVI    = self.InfectionRate*DVE[n] - self.RecoveryRate*DVI[n] - self.DeathRateCOVIDDV*DVI[n]

            #dSHR    = self.BirthRateHR*self.PopulationTotal - self.ContactRateSHR*(C/self.PopulationTotal)*SHR[n] + self.ImmunityLossRate*RHR[n] - self.VaccinationRateSHR*SHR[n]  - self.DeathRate*SHR[n] 
            #dEHR    = self.ContactRateSHR*(C/self.PopulationTotal)*SHR[n] - self.InfectionRate*EHR[n] - self.DeathRate*EHR[n] 
            #dIHR    = self.InfectionRate*EHR[n] - self.RecoveryRate*IHR[n] - self.DeathRateCOVIDHR*IHR[n]
            #dRHR    = self.RecoveryRate*IHR[n] - self.VaccinationRateSHR*RHR[n] - self.ImmunityLossRate*RHR[n] - self.DeathRate*RHR[n]

            #dSVHR     = self.VaccinationRateSHR*(SHR[n] + RHR[n]) - self.ContactRateSVHR*(C/self.PopulationTotal)*SVHR[n] - self.DeathRate*SVHR[n] - self.VaccinationRateSVHR*SVHR[n] + self.RecoveryRate*SVIHR[n]
            #dSVIHR    = self.InfectionRate*SVEHR[n] - self.RecoveryRate*SVIHR[n] - self.DeathRateCOVIDSVHR*SVIHR[n]
            #dSVEHR    = self.ContactRateSVHR*(C/self.PopulationTotal)*SVHR[n] - self.InfectionRate*SVEHR[n] - self.DeathRate*SVEHR[n]

            #dDVHR     = self.VaccinationRateSVHR*SVHR[n] - self.ContactRateDV*(C/self.PopulationTotal)*DVHR[n] + self.RecoveryRate*DVIHR[n] - self.DeathRate*DVHR[n] 
            #dDVEHR    = self.ContactRateDV*(C/self.PopulationTotal)*DVHR[n] - self.InfectionRate*DVEHR[n] - self.DeathRate*DVEHR[n]
            #dDVIHR    = self.InfectionRate*DVEHR[n] - self.RecoveryRate*DVIHR[n] - self.DeathRateCOVIDDVHR*DVIHR[n]

            S.append(S[n] + dS); E.append(E[n] + dE); I.append(I[n] + dI); R.append(R[n] + dR)
            SV.append(SV[n] + dSV); #SVE.append(SVE[n] + dSVE); SVI.append(SVI[n] + dSVI)
            DV.append(DV[n] + dDV);# DVE.append(SVE[n] + dDVE); DVI.append(SVI[n] + dDVI)
            #SHR.append(SHR[n] + dSHR); EHR.append(EHR[n] + dEHR); IHR.append(IHR[n] + dIHR); RHR.append(RHR[n] + dRHR)
            #SVHR.append(SVHR[n] + dSVHR); SVEHR.append(SVEHR[n] + dSVEHR); SVIHR.append(SVIHR[n] + dSVIHR)
            #DVHR.append(DVHR[n] + dDVHR); DVEHR.append(DVEHR[n] + dDVEHR); DVIHR.append(DVIHR[n] + dDVIHR)

            self.PopulationLR += dR + dS + dI + dE + dSV + dDV #+ dSVE + dSVI #+ dDVE + dDVI 
            #self.PopulationHR += dSHR + dEHR + dIHR + dRHR + dSVHR + dSVEHR + dSVIHR + dDVHR + dDVEHR + dDVIHR
            #self.PopulationTotal += self.PopulationLR #+ self.PopulationHR
        

        plt.plot(S, label = "Susceptible", c = '#1e68b3')
        plt.plot(np.array(I) + np.array(E), label = "Infected", c = '#b81111')
        plt.plot(R, label = "Recovered", c = '#aaacad')
        plt.plot(np.array(SV), label = "Vaccinated", c='#5e017d')
        plt.plot(np.array(DV), label = "Fully Vaccinated", c='#439603')

        #plt.plot(np.array(S) + np.array(SHR), label = "Susceptible", c='#1e68b3')
        #plt.plot(np.array(E) + np.array(EHR) + np.array(SVE) + \
        #         np.array(DVE) + np.array(SVEHR) + np.array(DVEHR) + \
        #         np.array(I) + np.array(IHR) + np.array(SVI) + \
        #         np.array(DVI) + np.array(DVIHR) + np.array(SVIHR), label="Infected", c='#b81111')
        #plt.plot(np.array(R) + np.array(RHR), label ="Recovered", c='#aaacad')
        #plt.plot(np.array(SV) + np.array(SVHR), label = "Vaccinated", c='#5e017d')
        #plt.plot(np.array(DV) + np.array(DVHR), label = "Fully Vaccinated", c='#439603')

        #plt.xlim(0,1000)
        #plt.ylim(0,1)
        plt.legend()
        plt.show()


if __name__ == "__main__":
    SEIRSV2 = SEIRSV2_Model()
    SEIRSV2.NumInt()