// Extract the tracks from a file using root command prompts
// example of use. The command
// >  root trackExtraction.C\(\"/Users/marcbergevin/RAT_ROOT/output.root\"\)
// Will apply this routine of the output.root file
// M.B
#include <iostream>
#include <iomanip>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClass.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>
#include <vector>
#include <TRandom3.h>

TRandom3 randNum;

//#include <libRATEvent.h>

int goldenFileExtractor(const char *file, const char *outfile = "null",
double nhit_min =3., double goodness_min = 0.1, double goodness_dir = 0.1,
double timeWindow_ns = 100000, double rate = 10.0, double maxDistance = 2.0,
double fidBoundR = 5.42, double fidBoundZ = 5.42,
double pmtBoundR = 6.42,double pmtBoundZ = 6.42,
double tankBoundR = 8.02635,double tankBoundZ = 8.02635) {

    // Define the incoming out outgoing Trees
    TFile *f = new TFile(file);
    TTree *tree = (TTree*) f->Get("T");
    if (tree==0x0){
        return -1;
    }
    TFile *f_out;
    printf("outfile: %s\n",outfile);
    if (TString(outfile) == TString("null")) {
        f_out = new TFile(Form("ntuple_%s",f->GetName()),"Recreate");
    }else{
        f_out = new TFile(outfile,"Recreate");
    }

    RAT::DS::Root *rds = new RAT::DS::Root();
    tree->SetBranchAddress("ds", &rds);
    int nEvents = tree->GetEntries();

    //Define all the analysis parameters
    Double_t totPE = 0.0,totNHIT = 0,totMom,goodness,dirGoodness,qTmp,timeTmp,timeTmp1,n9;
    Double_t oldX=-1e9,oldY=-1e9,oldZ=-1e9,newX,newY,newZ,dirX,dirY,dirZ,oldFVX=-1e9,oldFVY=-1e9,oldFVZ=-1e9,oldFVT;
    Double_t totQB = 0.0, q2 = 0.0, pmtCount = 0.0;
    Int_t ibd=0,es=0,cc=0,icc=0,nc=0,old_singal,evt;
    Double_t cosTheta,cosThetaSN,cosThetaSNIBD, local_time,local_time_tmp,delta_time,mc_nu_energy,mc_energy;

    int sub_ev=0,cnt_all = 0;
    Int_t subevents = 0,cnt,cntLoop;
    Int_t candidate,_candidate;
    Double_t timeDiff,r,z,r_t,z_t;

    TVector3 mcmomv_nu, mcmomv_particle;
    TVector3 posTruth,posReco,dirTruth,dirNu,dirReco,dirIBD,pos1,pos2;

    Int_t totcandidates, totMultiples;
    Double_t inner_dist,inner_time,inner_dist_fv,inner_time_fv;
    Int_t SV,old_FV,FV,GSV,IV,EV,OV,FV_t,GSV_t,IV_t,EV_t,OV_t,cnt_1,tot_FV,consecutive_FV;

    TTree *data = new TTree("data","low-energy detector events");
    data->Branch("pe",&totPE,"pe/D");
    data->Branch("nhit",&totNHIT,"nhit/D");
    data->Branch("n9",&n9,"n9/D");

    data->Branch("delta_time_ns",&delta_time,"delta_time_s/D");
    data->Branch("detected_ev",&sub_ev,"detected_ev/I");
    data->Branch("detected_ev_tot",&cnt,"detected_ev_tot/I");
    data->Branch("all_ev",&cnt_1,"all_ev/I");
    data->Branch("all_ev_tot",&cnt_all,"all_ev_tot/I");
    data->Branch("subevents",&subevents,"subevents/I");
    data->Branch("event_number",&evt,"event_number/I");

    data->Branch("candidate",&candidate,"candidate/I");

    data->Branch("mc_prim_energy",&mc_nu_energy,"mc_prim_energy/D");

    data->Branch("pos_goodness",&goodness,"pos_goodness/D");
    data->Branch("posReco","TVector3",&posReco,32000,0);
    data->Branch("reco_r",&r,"reco_r/D");
    data->Branch("reco_z",&z,"reco_z/D");

    data->Branch("posTruth","TVector3",&posTruth,32000,0);
    data->Branch("true_r",&r_t,"true_r/D");
    data->Branch("true_z",&z_t,"true_z/D");

    data->Branch("dir_goodness",&dirGoodness,"dir_goodness/D");
    data->Branch("dirReco","TVector3",&dirReco,32000,0);
    data->Branch("dirPrimaryMC","TVector3",&dirNu,32000,0);

    data->Branch("FV",&FV,"FV/I");
    data->Branch("GSV",&GSV,"GSV/I");
    data->Branch("EV",&EV,"EV/I");
    data->Branch("OV",&OV,"OV/I");
    data->Branch("IV",&IV,"IV/I");

    data->Branch("FV_truth",&FV_t,"FV_truth/I");
    data->Branch("GSV_truth",&GSV_t,"GSV_truth/I");
    data->Branch("EV_truth",&EV_t,"EV_truth/I");
    data->Branch("OV_truth",&OV_t,"OV_truth/I");
    data->Branch("IV_truth",&IV_t,"IV_truth/I");

    data->Branch("inner_dist",&inner_dist,"inner_dist/D");
    data->Branch("inner_time",&inner_time,"inner_time/D");
    data->Branch("inner_dist_fv",&inner_dist_fv,"inner_dist_fv/D");
    //    data->Branch("inner_time_fv",&inner_time_fv,"inner_time_fv/D");    //    data->Branch("old_FV",&old_FV,"old_FV/I");
    data->Branch("tot_FV",&tot_FV,"tot_FV/I");
    data->Branch("consecutive_FV",&consecutive_FV,"consecutive_FV/I");
    vector <double> subeventInfo;
    vector<vector <double> > eventInfo;
    RAT::DS::MC *mc;
    RAT::DS::MCParticle *prim;
    Int_t particleCountMC;
    for (evt = 0; evt < nEvents; evt++) {

        tree->GetEntry(evt);
        mc                          = rds->GetMC();
        particleCountMC             = mc->GetMCParticleCount();

        //Get the direction of the neutrino. Saved as last particle
        prim                        = mc->GetMCParticle(0);
        mcmomv_nu                   = prim->GetMomentum();
        dirNu                       = prim->GetMomentum();
        mc_nu_energy                = prim->ke;
        posTruth                    = prim->GetPosition();

        r_t = sqrt(pow(prim->GetPosition().X(),2)+ pow(prim->GetPosition().Y(),2))/1000.;
        z_t = prim->GetPosition().Z()/1000.;

        FindVolume(r_t,z_t,FV_t,GSV_t,IV_t,EV_t,OV_t,fidBoundR,fidBoundZ,pmtBoundR,pmtBoundZ,tankBoundR,tankBoundZ);

        //Find out how many subevents:
        subevents                   = rds->GetEVCount();
        timeTmp1                    = 0.;
        cnt = 0;

        subeventInfo.resize(0);
        eventInfo.resize(0);

        cnt_all = cnt = old_FV = tot_FV =  0 ;

        if(subevents == 0){
            //            resetVector(subeventInfo,eventInfo);
            sub_ev      = 0;
            inner_time  = 0;
            inner_dist  = 0;
            candidate   = 0;
            totNHIT     = 0;
            totPE       = 0;
            FV          = 0;
            GSV         = 0;
            IV          = 0;
            EV          = 0;
            OV          = 0;
            goodness    = 0;
            dirGoodness = 0;
            old_FV      = 0;
            r           = 0;
            z           = 0;
            posReco     = TVector3(0.,0.,0.);
            dirReco     = TVector3(0.,0.,0.);
            cnt_1       = 0;
            n9          = 0;
            data->Fill();

        }else{

            for (int k = 0; k<subevents; k++) {
                RAT::DS::EV *ev         = rds->GetEV(k);
                qTmp                    = ev->Nhits();
                RAT::DS::BonsaiFit *pb  = ev->GetBonsaiFit();
                goodness                = pb->GetGoodness();
                dirGoodness             = pb->GetDirGoodness();
                posReco                 = pb->GetPosition();
                totNHIT                 = pb->GetIDHit();
                totPE                   = pb->GetIDCharge();
                dirReco                 = pb->GetDirection();
                n9                      = pb->GetN9();
                newX                    = posReco.X();
                newY                    = posReco.Y();
                newZ                    = posReco.Z();
                dirX                    = dirReco.X();
                dirY                    = dirReco.Y();
                dirZ                    = dirReco.Z();

                timeTmp = timeTmp1 = 0;
                timeTmp = ev->GetCalibratedTriggerTime(); // 0 for first subevent, detlta for all others

                r = sqrt(pow(posReco.X(),2)+ pow(posReco.Y(),2))/1000.;
                z = posReco.Z()/1000.;

                if (n9 > nhit_min && goodness > goodness_min && dirGoodness > goodness_dir) {
                    cnt                  +=1;
                    cnt_all              +=1;

                    FindVolume(r,z,FV,GSV,IV,EV,OV,fidBoundR,fidBoundZ,pmtBoundR,pmtBoundZ,tankBoundR,tankBoundZ);
                    if (FV==1) {
                        tot_FV+=1;
                        consecutive_FV+=1;
                        inner_dist_fv =sqrt(pow(posReco.X()-oldFVX,2)+ pow(posReco.Y()-oldFVY,2)+pow(posReco.Z()-oldFVZ,2))/1000.;
                        oldFVX = posReco.X();
                        oldFVY = posReco.Y();
                        oldFVZ = posReco.Z();
                        //                    cout << inner_dist_fv << endl;


                    }else{
                        consecutive_FV=0;
                        inner_dist_fv = 0.;
                    }


                    //                printf("%d %f\n",cnt,timeTmp);
                    if (cnt ==1) {
                        inner_dist = sqrt(pow(newX-oldX,2)+ pow(newY-oldY,2) + pow(newZ-oldZ,2))/1000.;
                        if (timeTmp>0) {
                            inner_time = findNextTime(rate); // Find a random time
                        }else{
                            inner_time = -timeTmp; // To accomodate strange rat-pac backward time convention for Chains
                        }
                        oldX = newX;
                        oldY = newY;
                        oldZ = newZ;

                        timeTmp1        = timeTmp;
                        subeventInfo.push_back(cnt);                //0 sub_ev
                        subeventInfo.push_back(inner_time);         //1 time difference
                        subeventInfo.push_back(inner_dist);         //2 distance
                        if (timeDiff < timeWindow_ns && inner_dist < maxDistance) {
                            //                        eventInfo[cnt-2][3] = 1 ; // Issue with this if no event in middle
                            subeventInfo.push_back(1);              //3a is a candidate for accidentals
                        }else{
                            subeventInfo.push_back(0);              //3b is not a candidate
                        }                    subeventInfo.push_back(totNHIT);            //4 Record nhit
                        subeventInfo.push_back(totPE);              //5 record photoelectrons
                        subeventInfo.push_back(FV);                 //6 is in Fiducial Volume
                        subeventInfo.push_back(GSV);                //7 is in GammaShield volume
                        subeventInfo.push_back(IV);                 //8 is in InnerVolume
                        subeventInfo.push_back(EV);                 //9 is in EntireVolume
                        subeventInfo.push_back(OV);                 //10 is in outerVolume
                        subeventInfo.push_back(goodness);           //11 bonsai position goodness
                        subeventInfo.push_back(dirGoodness);        //12 bonsai direction goodness
                        subeventInfo.push_back(old_FV);             //13 Previous event was in FV
                        subeventInfo.push_back(r);                  //14 Reconstructed radius
                        subeventInfo.push_back(z);                  //15 reconstruced z position
                        subeventInfo.push_back(newX);               //16 New x coordinate
                        subeventInfo.push_back(newY);               //17 New y coordinate
                        subeventInfo.push_back(newZ);               //18 New z coordinate
                        subeventInfo.push_back(dirX);               //19 New x direction
                        subeventInfo.push_back(dirY);               //20 New y direction
                        subeventInfo.push_back(dirZ);               //21 New z direction
                        subeventInfo.push_back(cnt_all);            //22 All sub-ev
                        subeventInfo.push_back(n9);                 //23 hit in 9 nanosecond

                        eventInfo.push_back(subeventInfo);

                        subeventInfo.resize(0);

                    }else if(cnt>1){
                        inner_dist = sqrt(pow(newX-oldX,2)+ pow(newY-oldY,2) + pow(newZ-oldZ,2))/1000.;
                        oldX = newX;
                        oldY = newY;
                        oldZ = newZ;

                        timeDiff   = timeTmp-timeTmp1;
                        timeTmp1        = timeTmp;
                        subeventInfo.push_back(cnt);                //0 sub_ev
                        subeventInfo.push_back(abs(timeDiff));      //1 time difference
                        subeventInfo.push_back(inner_dist);         //2 distance
                        if (timeDiff < timeWindow_ns && inner_dist < maxDistance) {
                            //                        eventInfo[cnt-2][3] = 1 ; // Issue with this if no event in middle
                            subeventInfo.push_back(1);              //3a is a candidate
                            if(eventInfo[cnt-2][13] == 1){
                                old_FV = 1;
                            }else{
                                old_FV = 0;
                            }
                        }else{
                            subeventInfo.push_back(0);              //3b is not a candidate
                        }
                        subeventInfo.push_back(totNHIT);            //4 Record nhit
                        subeventInfo.push_back(totPE);              //5 record photoelectrons
                        subeventInfo.push_back(FV);                 //6 is in Fiducial Volume
                        subeventInfo.push_back(GSV);                //7 is in GammaShield volume
                        subeventInfo.push_back(IV);                 //8 is in InnerVolume
                        subeventInfo.push_back(EV);                 //9 is in EntireVolume
                        subeventInfo.push_back(OV);                 //10 is in outerVolume
                        subeventInfo.push_back(goodness);           //11 bonsai position goodness
                        subeventInfo.push_back(dirGoodness);        //12 bonsai direction goodness
                        subeventInfo.push_back(old_FV);             //13 Previous event was in FV
                        subeventInfo.push_back(r);                  //14 Reconstructed radius
                        subeventInfo.push_back(z);                  //15 reconstruced z position
                        subeventInfo.push_back(newX);               //16 New x coordinate
                        subeventInfo.push_back(newY);               //17 New y coordinate
                        subeventInfo.push_back(newZ);               //18 New z coordinate
                        subeventInfo.push_back(dirX);               //19 New x direction
                        subeventInfo.push_back(dirY);               //20 New y direction
                        subeventInfo.push_back(dirZ);               //21 New z direction
                        subeventInfo.push_back(cnt_all);            //22 All sub-ev
                        subeventInfo.push_back(n9);                 //23 hit in 9 nanosecond

                        eventInfo.push_back(subeventInfo);
                        subeventInfo.resize(0);

                    }
                }else{
                    cnt_all     +=1;
                    subeventInfo.push_back(0);                      //0 sub_ev
                    subeventInfo.push_back(0.0);                   //1 time difference
                    subeventInfo.push_back(0.0);                   //2 distance
                    subeventInfo.push_back(0);                      //3 is a candidate
                    subeventInfo.push_back(totNHIT);                //4 Record nhit
                    subeventInfo.push_back(totPE);                  //5 record photoelectrons
                    subeventInfo.push_back(0);                     //6 is in Fiducial Volume

                    subeventInfo.push_back(0);                      //7 is in GammaShield volume
                    subeventInfo.push_back(0);                     //8 is in InnerVolume
                    subeventInfo.push_back(0);                     //9 is in EntireVolume
                    subeventInfo.push_back(0);                     //10 is in outerVolume
                    subeventInfo.push_back(goodness);               //11 bonsai position goodness
                    subeventInfo.push_back(dirGoodness);            //12 bonsai direction goodness
                    subeventInfo.push_back(0);                      //13 Previous event was in FV
                    subeventInfo.push_back(r);                      //14 Reconstructed radius
                    subeventInfo.push_back(z);                      //15 reconstruced z position
                    subeventInfo.push_back(newX);                   //16 New x coordinate
                    subeventInfo.push_back(newY);                   //17 New y coordinate
                    subeventInfo.push_back(newZ);                   //18 New z coordinate
                    subeventInfo.push_back(dirX);                   //19 New x direction
                    subeventInfo.push_back(dirY);                   //20 New y direction
                    subeventInfo.push_back(dirZ);                   //21 New z direction
                    subeventInfo.push_back(cnt_all);              //22 All sub-ev
                    subeventInfo.push_back(n9);                 //23 hit in 9 nanosecond

                    eventInfo.push_back(subeventInfo);
                    subeventInfo.resize(0);

                }
            }

        }

        for (unsigned long evt_idx = 0; evt_idx < eventInfo.size();evt_idx++) {
            sub_ev      = eventInfo[evt_idx][0];
            inner_time  = eventInfo[evt_idx][1];
            inner_dist  = eventInfo[evt_idx][2];
            candidate   = eventInfo[evt_idx][3];
            totNHIT     = eventInfo[evt_idx][4];
            totPE       = eventInfo[evt_idx][5];
            FV          = eventInfo[evt_idx][6];
            GSV         = eventInfo[evt_idx][7];
            IV          = eventInfo[evt_idx][8];
            EV          = eventInfo[evt_idx][9];
            OV          = eventInfo[evt_idx][10];
            goodness    = eventInfo[evt_idx][11];
            dirGoodness = eventInfo[evt_idx][12];
            old_FV      = eventInfo[evt_idx][13];
            r           = eventInfo[evt_idx][14];
            z           = eventInfo[evt_idx][15];
            posReco     = TVector3(eventInfo[evt_idx][16],eventInfo[evt_idx][17],eventInfo[evt_idx][18]);
            dirReco     = TVector3(eventInfo[evt_idx][19],eventInfo[evt_idx][20],eventInfo[evt_idx][21]);
            cnt_1       = eventInfo[evt_idx][22];
            n9          = eventInfo[evt_idx][23];
            data->Fill();
        }

    }

    data->Write();
    rds->Delete();
    f_out->Close();
    f->Close();
    return 0;
}

void FindVolume(Double_t r_t,Double_t z_t,Int_t &FV_t,Int_t &GSV_t,Int_t &IV_t,Int_t &EV_t,Int_t &OV_t,double fidBoundR,double fidBoundZ,double pmtBoundR,double pmtBoundZ, double tankBoundR, double tankBoundZ){

    FV_t = GSV_t = IV_t = EV_t = OV_t = 0;
    if (r_t < fidBoundR && fabs(z_t) < fidBoundZ) {
        FV_t =1;
        IV_t =1;
        EV_t =1;
    }else if(!(r_t < fidBoundR && fabs(z_t) < fidBoundZ) && (r_t < pmtBoundR && fabs(z_t) < pmtBoundZ)){
        GSV_t =1;
        IV_t  =1;
        EV_t =1;
    }else if(!(r_t < pmtBoundR && fabs(z_t) < pmtBoundZ) && (r_t < tankBoundR && fabs(z_t) < tankBoundZ)){
        OV_t =1;
        EV_t =1;
    }else{
        EV_t = 0;
    }
}

Double_t findNextTime(Double_t rate){

    return randNum.Exp(1./rate)*1.0e9;
}

void resetVector(    vector <double> subeventInfo,vector<vector <double> > &eventInfo){
    eventInfo.resize(0);
    subeventInfo.resize(0);
    subeventInfo.push_back(0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    subeventInfo.push_back(0.0);
    eventInfo.push_back(subeventInfo);


}
