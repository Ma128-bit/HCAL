import ROOT
import time
import sys
import os
c_func = """
double langaufun(double *x, double *par) {
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location
 
      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
 
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
 
      return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1F *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi,
double *fitparams, double *fiterrors, double *ChiSqr, int *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf
 
   int i;
   char FunName[100];
 
   snprintf(FunName, sizeof(FunName), "Fitfcn_%s", his->GetName());
 
   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;
 
   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");
 
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }
 
   his->Fit(FunName,"QRB0");   // fit within specified range, use ParLimits, do not plot
 
   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf
 
   return (ffit);              // return fit function
 
}

int langaupro(double *params, double &maxx, double &FWHM) {
 
   // Searches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.
 
   double p,x,fy,fxr,fxl;
   double step;
   double l,lold;
   int i = 0;
   int MAXCALLS = 10000;
 
 
   // Search for maximum
 
   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;
 
 
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
 
      lold = l;
      x = p + step;
      l = langaufun(&x,params);
 
      if (l < lold)
         step = -step/10;
 
      p += step;
   }
 
   if (i == MAXCALLS)
      return (-1);
 
   maxx = x;
 
   fy = l/2;
 
 
   // Search for right x location of fy
 
   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;
 
 
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
 
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }
 
   if (i == MAXCALLS)
      return (-2);
 
   fxr = x;
 
 
   // Search for left x location of fy
 
   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;
 
   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
 
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }
 
   if (i == MAXCALLS)
      return (-3);
 
 
   fxl = x;
 
   FWHM = fxr - fxl;
   return (0);
}

vector<vector<short>> new_raw_q(unsigned long long evt, vector<vector<short>> raw_q, vector<string> mmChamber, vector<int> mmStrip, vector<int> t_max_q){
    vector<vector<short>> raw_q_temp = raw_q;
    cout<<"Event: "<<evt<<endl;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        if(mmChamber.at(i) == "Tmm-Lecce2" || mmChamber.at(i) == "Tmm-Lecce1"){
            TString st=Form("event_%llu_strip_%d_onChamber_",evt, mmStrip.at(i));
            TH1F *q_per_strip = new TH1F(st+(TString)mmChamber.at(i), st+(TString)mmChamber.at(i), 15, 0, 15);
            double fr[2];

            Int_t j=0;
            bool sel_mim=false;
            double charge_max= 0.0;
            double qtot=0.0;
            for(j=0; j<raw_q.at(i).size() ; j++){
                double charge= (raw_q.at(i)).at(j);
                qtot=qtot+charge;
                if(charge_max<charge){
                    charge_max=charge;
                }
                if(sel_mim==false & charge>1){
                    fr[0]=j+1;
                    sel_mim=true;
                }
                if(charge>1){
                    fr[1]=j;
                }
                if (charge<=0){
                    charge = 0.001;
                }
                q_per_strip->SetBinContent(j+1, charge);
            }

            double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
                
            sv[0]=TMath::Sqrt(q_per_strip->GetRMS()); sv[1]=q_per_strip->GetMean(); sv[2]=qtot; sv[3]=TMath::Sqrt(q_per_strip->GetRMS());
            pllo[0]=sv[0]/4; pllo[1]=sv[1]/4; pllo[2]=sv[2]/4; pllo[3]=sv[3]/4;
            plhi[0]=sv[0]*4; plhi[1]=sv[1]*4; plhi[2]=sv[2]*1000; plhi[3]=sv[3]*4;
                
            double chisqr;
            int    ndf;
         
            if(charge_max>40){
                TF1 *fitsnr = langaufit(q_per_strip,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
            
                delete q_per_strip;
            
                Double_t chiSquare = chisqr;
                int Ndof_ToT = ndf;
            
                if(chiSquare/(Ndof_ToT)<6 && Ndof_ToT>1 && fp[0]<5.5 && fp[3]<5.5){
                        //cout<<"Larhezza landau: "<<fp[0]<<" MPV: "<<fp[1]<<"--"<<" Larhezza gauss: "<<fp[3]<<endl;
                        //c->cd();
                        //q_per_strip->Draw();
                        //fitsnr->Draw("same");
                        //c->SaveAs("Plot/"+st+(TString)mmChamber.at(i)+".png");
                        //cout<<st+(TString)mmChamber.at(i)<<endl;
                        //cout<<"chiSquare: "<<chiSquare<<endl;
                        //cout<<"ndof: "<<Ndof_ToT<<endl;
                        //c->Clear();
                }
                
                else{
                    for(j=0; j<raw_q.at(i).size() ; j++){
                        (raw_q_temp.at(i)).at(j) = 0.0;
                    }
                }
                delete fitsnr;
            }
            else{
                delete q_per_strip;
                for(j=0; j<raw_q.at(i).size() ; j++){
                    (raw_q_temp.at(i)).at(j) = 0.0;
                }
            }
        }
            
    }//End mmStrip.size() loop
    return raw_q_temp;
}


vector<vector<short>> histo_var(unsigned long long evt, vector<vector<short>> raw_q, vector<string> mmChamber, vector<int> mmStrip, vector<int> t_max_q){
    vector<vector<short>> raw_q_temp = raw_q;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        vector<double> histo_temp = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        if(mmChamber.at(i) == "Tmm-Lecce2" || mmChamber.at(i) == "Tmm-Lecce1"){
            TString st=Form("event_%llu_strip_%d_onChamber_",evt, mmStrip.at(i));
            TH1F *q_per_strip = new TH1F(st+(TString)mmChamber.at(i), st+(TString)mmChamber.at(i), 15, 0, 15);
            double meanY=0.0;
            int count=0;
            for(Int_t j=0; j<raw_q.at(i).size() ; j++){
                double charge= (raw_q.at(i)).at(j);
                meanY = meanY+charge;
                count++;
                //if (charge<=0){ charge = 0.001; }
                q_per_strip->SetBinContent(j+1, charge);
            }
            meanY=meanY/count;
            double SquaredDevY = 0.0;
            for(Int_t j=0; j<raw_q.at(i).size() ; j++){
                double charge= (raw_q.at(i)).at(j);
                double deviation = charge - meanY;
                SquaredDevY += deviation * deviation;
            }
            SquaredDevY=TMath::Sqrt(SquaredDevY / count);
            histo_temp[0]=q_per_strip->GetMean();
            histo_temp[1]=meanY;
            histo_temp[2]=TMath::Sqrt(q_per_strip->GetRMS());
            histo_temp[3]=SquaredDevY;
            histo_temp[4]=q_per_strip->GetKurtosis();
            histo_temp[5]=q_per_strip->GetSkewness();
            histo_temp[6]=q_per_strip->GetMaximumBin();
            if(histo_temp[0]<2 || histo_temp[2]<1 || histo_temp[1]<2 || histo_temp[5]>1.5 || histo_temp[4]>3 || histo_temp[3]<2 || (histo_temp[1]>400 && histo_temp[3]<100)){
                for(Int_t j=0; j<raw_q.at(i).size() ; j++){
                    (raw_q_temp.at(i)).at(j) = 0.0;
                }
            }
            delete q_per_strip;
        }
    }//End mmStrip.size() loop
    return raw_q_temp;
}

vector<short> new_max_q(vector<short> max_q, vector<vector<short>> raw_q, vector<string> mmChamber, vector<int> mmStrip, vector<int> t_max_q){
    vector<short> max_q_temp = max_q;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        if(mmChamber.at(i) == "Tmm-Lecce2" || mmChamber.at(i) == "Tmm-Lecce1"){
            short max_charge = 0.0;
            for(Int_t j=0; j<raw_q.at(i).size() ; j++){
                short charge= (raw_q.at(i)).at(j);
                if(charge > max_charge){
                    max_charge = charge;
                }
            }
            max_q_temp.at(i) = max_charge;
        }
    }//End mmStrip.size() loop
    return max_q_temp;
}

vector<int> new_t_max_q(vector<short> max_q, vector<vector<short>> raw_q, vector<string> mmChamber, vector<int> mmStrip, vector<int> t_max_q){
    vector<int> t_max_q_temp = t_max_q;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        if(mmChamber.at(i) == "Tmm-Lecce2" || mmChamber.at(i) == "Tmm-Lecce1"){
            Int_t j=0;
            short max_charge = 0.0;
            int temp_max = 0.0;
            for(j=0; j<raw_q.at(i).size() ; j++){
                short charge= (raw_q.at(i)).at(j);
                if(charge > max_charge){
                    max_charge = charge;
                    temp_max = j+1;
                }
            }
            t_max_q_temp.at(i) = temp_max;
        }
    }//End mmStrip.size() loop
    return t_max_q_temp;
}

vector<vector<short>> conv_raw_q(ROOT::VecOps::RVec<vector<short>> raw_q, ROOT::VecOps::RVec<int> mmStrip){
    vector<vector<short>> temp;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        vector<short> pippo;
        for(Int_t j=0; j<raw_q.at(i).size() ; j++){
            pippo.push_back((short)(raw_q.at(i).at(j)));
        }
        temp.push_back(pippo);
    }
    return temp;
}

vector<string> conv_mmChamber(ROOT::VecOps::RVec<string> mmChamber, ROOT::VecOps::RVec<int> mmStrip){
    vector<string> temp;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        temp.push_back(mmChamber.at(i));
    }
    return temp;
}

vector<unsigned int> conv_uint(ROOT::VecOps::RVec<unsigned int> x_var, ROOT::VecOps::RVec<int> mmStrip){
    vector<unsigned int> temp;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        temp.push_back(x_var.at(i));
    }
    return temp;
}

vector<int> conv_int(ROOT::VecOps::RVec<int> x_var, ROOT::VecOps::RVec<int> mmStrip){
    vector<int> temp;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        temp.push_back(x_var.at(i));
    }
    return temp;
}

vector<char> conv_char(ROOT::VecOps::RVec<char> x_var, ROOT::VecOps::RVec<int> mmStrip){
    vector<char> temp;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        temp.push_back(x_var.at(i));
    }
    return temp;
}

vector<short> conv_short(ROOT::VecOps::RVec<short> x_var, ROOT::VecOps::RVec<int> mmStrip){
    vector<short> temp;
    for(Int_t i = 0; i<mmStrip.size() ; i++){
        temp.push_back(x_var.at(i));
    }
    return temp;
}

"""
start_time = time.time() #Misura tempo

ROOT.gInterpreter.Declare(c_func) # Dichiarare la funzione all'interprete ROOT

#ROOT.EnableImplicitMT() #Per eseguire in parallelo

args = sys.argv
if len(args) != 3:
    print("the argument should be ./run#(in) and out")
    print("Length: ",len(args))
name = args[1]
name2 = args[2]

current_directory = os.getcwd()
dictionaries_path = "./dictionaries"
os.chdir(dictionaries_path)
ROOT.gInterpreter.GenerateDictionary("vector<vector<short> >","vector")
ROOT.gInterpreter.GenerateDictionary("vector<string>","vector")
print("End Generation of Dictionaries")

df = ROOT.RDataFrame("apv_raw", name+".root")

#Conversione variabili:
df = df.Redefine("raw_q","conv_raw_q(raw_q,mmStrip)")
df = df.Redefine("mmChamber","conv_mmChamber(mmChamber,mmStrip)")
df = df.Redefine("srsFec","conv_uint(srsFec,mmStrip)")
df = df.Redefine("srsFec","conv_uint(srsFec,mmStrip)")
df = df.Redefine("srsChip","conv_uint(srsChip,mmStrip)")
df = df.Redefine("srsChan","conv_uint(srsChan,mmStrip)")
df = df.Redefine("mmLayer","conv_int(mmLayer,mmStrip)")
df = df.Redefine("t_max_q","conv_int(t_max_q,mmStrip)")
df = df.Redefine("mmReadout","conv_char(mmReadout,mmStrip)")
df = df.Redefine("max_q","conv_short(max_q,mmStrip)")
df = df.Redefine("mmStrip","conv_int(mmStrip,mmStrip)")

df = df.Redefine("raw_q", "histo_var(evt, raw_q, mmChamber, mmStrip, t_max_q)")
df = df.Redefine("max_q", "new_max_q(max_q, raw_q, mmChamber, mmStrip, t_max_q)")
df = df.Redefine("t_max_q", "new_t_max_q(max_q, raw_q, mmChamber, mmStrip, t_max_q)")

df = df.Redefine("raw_q","new_raw_q(evt, raw_q, mmChamber, mmStrip, t_max_q)")
df = df.Redefine("max_q", "new_max_q(max_q, raw_q, mmChamber, mmStrip, t_max_q)")
df = df.Redefine("t_max_q", "new_t_max_q(max_q, raw_q, mmChamber, mmStrip, t_max_q)")

new_tree = df.Snapshot("apv_raw", name2+"_new2.root") # Creare un nuovo file ROOT con un nuovo TTree

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Tempo trascorso: {elapsed_time} secondi")
