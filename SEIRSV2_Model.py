import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class SEIRSV2_Model():

    def __init__(self):

        COVIDdf = (pd.read_csv('CSVfiles/owid-covid-data-uk.csv')).replace(np.nan, 0)

        #General Coefficients
        self.Population             = 1
        self.InitialExposures       = 0.0000056
        self.PopulationHR           = 0
        self.VaccineDay             = 0
        self.BirthRate              = 2.78363e-5 
        self.DeathRate              = 2.48321e-5
        
        #Pandemic Coefficients
        self.InfectionRate          = 0.19
        self.RecoveryRate           = 1/10
        self.ImmunityLossRate       = 1/210
        self.ContactRateS           = (self.RecoveryRate)*np.array(COVIDdf['reproduction_rate'][30:1030])
        self.ContactRateSV          =  self.ContactRateS * 0.19
        self.ContactRateDV          =  self.ContactRateS * 0.09
        self.ContactRateSHR         = (self.RecoveryRate)*np.array(COVIDdf['reproduction_rate'][30:1030]) * 2
        self.ContactRateSVHR        =  self.ContactRateSHR * 0.19
        self.ContactRateDVHR        =  self.ContactRateSHR * 0.09
        self.VaccinationRateS       = 0
        self.VaccinationRateSV      = 0
        self.VaccinationRateSHR     = 0
        self.VaccinationRateSVHR    = 0
        self.DeathRateCOVID         = 8.95027e-4
        self.DeathRateCOVIDHR       = self.DeathRateCOVID * 2
        self.DeathRateCOVIDSV       = self.DeathRateCOVID * 0.19
        self.DeathRateCOVIDDV       = self.DeathRateCOVID * 0.09
        self.DeathRateCOVIDSVHR     = self.DeathRateCOVIDHR * 0.19
        self.DeathRateCOVIDDVHR     = self.DeathRateCOVIDHR * 0.09

    def NumInt(self):

        S = [self.Population - self.InitialExposures]; E = [self.InitialExposures]; I = [0]; R = [0]
        SV  = [0]; SVE = [0]; SVI = [0]
        DV  = [0]; DVE = [0]; DVI = [0]
        SHR = [self.PopulationHR]; EHR = [0]; IHR = [0]; RHR = [0]
        SVHR = [0]; SVEHR = [0]; SVIHR = [0]
        DVHR = [0]; DVEHR = [0]; DVIHR = [0]

        for n in np.arange(0, 1000):

            C = I[n] + E[n] + SVI[n] + SVE[n] + DVI[n] + DVE[n] + EHR[n] + IHR[n] + SVIHR[n] + SVEHR[n] + DVEHR[n] + DVIHR[n] 

            dS      = self.BirthRate*self.Population - self.ContactRateS[n]*(C/self.Population)*S[n] + self.ImmunityLossRate*R[n] - self.VaccinationRateS*S[n]  - self.DeathRate*S[n] 
            dE      = self.ContactRateS[n]*(C/self.Population)*S[n] - self.InfectionRate*E[n] - self.DeathRate*E[n] 
            dI      = self.InfectionRate*E[n] - self.RecoveryRate*I[n] - self.DeathRateCOVID*I[n]
            dR      = self.RecoveryRate*I[n] - self.VaccinationRateS*R[n] - self.ImmunityLossRate*R[n] - self.DeathRate*R[n]

            dSV     = self.VaccinationRateS*(S[n] + R[n]) - self.ContactRateSV[n]*(C/self.Population)*SV[n] - self.DeathRate*SV[n] - self.VaccinationRateSV*SV[n] + self.RecoveryRate*SVI[n]
            dSVE    = self.ContactRateSV[n]*(C/self.Population)*SV[n] - self.InfectionRate*SVE[n] - self.DeathRate*SVE[n]
            dSVI    = self.InfectionRate*SVE[n] - self.RecoveryRate*SVI[n] - self.DeathRateCOVIDSV*SVI[n]

            dDV     = self.VaccinationRateSV*SV[n] - self.ContactRateDV[n]*(C/self.Population)*DV[n] + self.RecoveryRate*DVI[n] - self.DeathRate*DV[n] 
            dDVE    = self.ContactRateDV[n]*(C/self.Population)*DV[n] - self.InfectionRate*DVE[n] - self.DeathRate*DVE[n]
            dDVI    = self.InfectionRate*DVE[n] - self.RecoveryRate*DVI[n] - self.DeathRateCOVIDDV*DVI[n]

            dSHR    = - self.ContactRateSHR[n]*(C/self.Population)*SHR[n] + self.ImmunityLossRate*RHR[n] - self.VaccinationRateSHR*SHR[n]  - self.DeathRate*SHR[n] 
            dEHR    = self.ContactRateSHR[n]*(C/self.Population)*SHR[n] - self.InfectionRate*EHR[n] - self.DeathRate*EHR[n] 
            dIHR    = self.InfectionRate*EHR[n] - self.RecoveryRate*IHR[n] - self.DeathRateCOVIDHR*IHR[n]
            dRHR    = self.RecoveryRate*IHR[n] - self.VaccinationRateSHR*RHR[n] - self.ImmunityLossRate*RHR[n] - self.DeathRate*RHR[n]

            dSVHR     = self.VaccinationRateSHR*(SHR[n] + RHR[n]) - self.ContactRateSVHR[n]*(C/self.Population)*SVHR[n] - self.DeathRate*SVHR[n] - self.VaccinationRateSVHR*SVHR[n] + self.RecoveryRate*SVIHR[n]
            dSVEHR    = self.ContactRateSVHR[n]*(C/self.Population)*SVHR[n] - self.InfectionRate*SVEHR[n] - self.DeathRate*SVEHR[n]
            dSVIHR    = self.InfectionRate*SVEHR[n] - self.RecoveryRate*SVIHR[n] - self.DeathRateCOVIDSVHR*SVIHR[n]

            dDVHR     = self.VaccinationRateSVHR*SVHR[n] - self.ContactRateDV[n]*(C/self.Population)*DVHR[n] + self.RecoveryRate*DVIHR[n] - self.DeathRate*DVHR[n] 
            dDVEHR    = self.ContactRateDV[n]*(C/self.Population)*DVHR[n] - self.InfectionRate*DVEHR[n] - self.DeathRate*DVEHR[n]
            dDVIHR    = self.InfectionRate*DVEHR[n] - self.RecoveryRate*DVIHR[n] - self.DeathRateCOVIDDVHR*DVIHR[n]

            S.append(S[n] + dS); E.append(E[n] + dE); I.append(I[n] + dI); R.append(R[n] + dR)
            SV.append(SV[n] + dSV); SVE.append(SVE[n] + dSVE); SVI.append(SVI[n] + dSVI)
            DV.append(SV[n] + dDV); DVE.append(SVE[n] + dDVE); DVI.append(SVI[n] + dDVI)
            SHR.append(SHR[n] + dSHR); EHR.append(EHR[n] + dEHR); IHR.append(IHR[n] + dIHR); RHR.append(RHR[n] + dRHR)
            SVHR.append(SVHR[n] + dSVHR); SVEHR.append(SVEHR[n] + dSVEHR); SVIHR.append(SVIHR[n] + dSVIHR)
            DVHR.append(DVHR[n] + dDVHR); DVEHR.append(DVEHR[n] + dDVEHR); DVIHR.append(DVIHR[n] + dDVIHR)

            self.Population += dR + dS + dI + dE + dSV + dSVE + dSVI + \
                               dDV + dDVE + dDVI + dSHR + dEHR + dIHR + \
                               dRHR + dSVHR + dSVEHR + dSVIHR + dDVHR + \
                               dDVEHR + dDVIHR
        
        plt.plot(I)
        plt.show()


if __name__ == "__main__":
    SEIRSV2 = SEIRSV2_Model()
    SEIRSV2.NumInt()

