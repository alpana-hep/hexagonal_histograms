#include "include/auxiliary_boundary_lines.h"
#include <map>

void beautify_plot(bool drawLine = true, bool drawText = true);

//    std::map<int, string> map_HGCROC_pin_str = {18:"CALIB0", 57:"CALIB1", 96:"CALIB0", 135:"CALIB1", 174:"CALIB0", 213:"CALIB1"}
// key=globalId, value=HGCROC_pin, "CALIB" is replaced with the same pin number as the companion cell
std::map<int, int> map_HGCROC_pin = {{0,0}, {1,1}, {2,2}, {3,3}, {4,4}, {5,5}, {6,6}, {7,7}, {9,9}, {10,10}, {11,11}, {12,12}, {13,13}, {14,14}, {15,15}, {16,16}, {20,19}, {21,20}, {22,21}, {23,22}, {24,23}, {25,24}, {26,25}, {27,26}, {29,28}, {30,29}, {31,30}, {32,31}, {33,32}, {34,33}, {35,34}, {36,35}, {39,36}, {40,37}, {41,38}, {42,39}, {43,40}, {44,41}, {45,42}, {46,43}, {48,45}, {49,46}, {50,47}, {51,48}, {52,49}, {53,50}, {54,51}, {55,52}, {59,55}, {60,56}, {61,57}, {62,58}, {63,59}, {64,60}, {65,61}, {66,62}, {68,64}, {69,65}, {70,66}, {71,67}, {72,68}, {73,69}, {74,70}, {75,71}, {78,0}, {79,1}, {80,2}, {81,3}, {82,4}, {83,5}, {84,6}, {85,7}, {87,9}, {88,10}, {89,11}, {90,12}, {91,13}, {92,14}, {93,15}, {94,16}, {98,19}, {99,20}, {100,21}, {101,22}, {102,23}, {103,24}, {104,25}, {105,26}, {107,28}, {108,29}, {109,30}, {110,31}, {111,32}, {112,33}, {113,34}, {114,35}, {117,36}, {118,37}, {119,38}, {120,39}, {121,40}, {122,41}, {123,42}, {124,43}, {126,45}, {127,46}, {128,47}, {129,48}, {130,49}, {131,50}, {132,51}, {133,52}, {137,55}, {138,56}, {139,57}, {140,58}, {141,59}, {142,60}, {143,61}, {144,62}, {146,64}, {147,65}, {148,66}, {149,67}, {150,68}, {151,69}, {152,70}, {153,71}, {156,0}, {157,1}, {158,2}, {159,3}, {160,4}, {161,5}, {162,6}, {163,7}, {165,9}, {166,10}, {167,11}, {168,12}, {169,13}, {170,14}, {171,15}, {172,16}, {176,19}, {177,20}, {178,21}, {179,22}, {180,23}, {181,24}, {182,25}, {183,26}, {185,28}, {186,29}, {187,30}, {188,31}, {189,32}, {190,33}, {191,34}, {192,35}, {195,36}, {196,37}, {197,38}, {198,39}, {199,40}, {200,41}, {201,42}, {202,43}, {204,45}, {205,46}, {206,47}, {207,48}, {208,49}, {209,50}, {210,51}, {211,52}, {215,55}, {216,56}, {217,57}, {218,58}, {219,59}, {220,60}, {221,61}, {222,62}, {224,64}, {225,65}, {226,66}, {227,67}, {228,68}, {229,69}, {230,70}, {231,71}, {18,12}, {57,55}, {96,12}, {135,55}, {174,12}, {213,55} };

// key=globalId, value=padId
std::map<int, int> map_SiCell_padId = {{0,36}, {1,26}, {2,35}, {3,25}, {4,8}, {5,17}, {6,16}, {7,7}, {9,6}, {10,15}, {11,5}, {12,13}, {13,34}, {14,33}, {15,23}, {16,24}, {18,14}, {20,4}, {21,12}, {22,3}, {23,11}, {24,2}, {25,1}, {26,10}, {27,9}, {29,21}, {30,31}, {31,22}, {32,32}, {33,19}, {34,29}, {35,20}, {36,30}, {39,46}, {40,58}, {41,47}, {42,59}, {43,44}, {44,57}, {45,56}, {46,45}, {48,74}, {49,88}, {50,87}, {51,73}, {52,71}, {53,86}, {54,72}, {55,85}, {57,70}, {59,69}, {60,84}, {61,68}, {62,83}, {63,43}, {64,42}, {65,55}, {66,54}, {68,67}, {69,82}, {70,66}, {71,81}, {72,52}, {73,41}, {74,53}, {75,40}, {78,165}, {79,177}, {80,151}, {81,166}, {82,178}, {83,188}, {84,198}, {85,189}, {87,179}, {88,167}, {89,168}, {90,153}, {91,137}, {92,123}, {93,138}, {94,152}, {96,154}, {98,155}, {99,139}, {100,140}, {101,125}, {102,126}, {103,111}, {104,110}, {105,95}, {107,109}, {108,93}, {109,124}, {110,108}, {111,94}, {112,79}, {113,80}, {114,65}, {117,121}, {118,135}, {119,136}, {120,150}, {121,107}, {122,106}, {123,91}, {124,122}, {126,120}, {127,104}, {128,89}, {129,105}, {130,76}, {131,75}, {132,90}, {133,60}, {135,62}, {137,61}, {138,48}, {139,49}, {140,37}, {141,92}, {142,78}, {143,77}, {144,63}, {146,38}, {147,27}, {148,28}, {149,18}, {150,39}, {151,64}, {152,50}, {153,51}, {156,99}, {157,98}, {158,115}, {159,114}, {160,96}, {161,97}, {162,113}, {163,112}, {165,127}, {166,128}, {167,141}, {168,142}, {169,130}, {170,145}, {171,144}, {172,129}, {174,143}, {176,157}, {177,156}, {178,169}, {179,170}, {180,180}, {181,190}, {182,181}, {183,191}, {185,171}, {186,172}, {187,158}, {188,159}, {189,192}, {190,193}, {191,182}, {192,183}, {195,100}, {196,101}, {197,116}, {198,117}, {199,146}, {200,132}, {201,147}, {202,131}, {204,102}, {205,103}, {206,119}, {207,118}, {208,148}, {209,134}, {210,133}, {211,149}, {213,163}, {215,162}, {216,164}, {217,175}, {218,176}, {219,160}, {220,173}, {221,161}, {222,174}, {224,186}, {225,187}, {226,196}, {227,197}, {228,195}, {229,184}, {230,185}, {231,194}};

void th2poly(TString inputfile,TString inputfile1, TString input_calib, TString png_String, TString outputfile, double range, bool drawLine=false, TString layer_="-1", TString runNumber="-1", int ped_method=-1)
{
    TCanvas *c1 = new TCanvas("c1", "", 900, 900);
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.1);
    if(ped_method==0 || ped_method==1 || ped_method==6)
      gStyle->SetPaintTextFormat(".0f");
    else
      gStyle->SetPaintTextFormat(".2f");
    // ----------------- //
    //reading the cell area and global channel ID //
    std::ifstream area_file;
    area_file.open("table_cell_area_in_LD_full_wafer_v1.txt",ios::in);
    std::vector<int> area_globalID;
    std::vector<int> area_cellArea;
    if(!area_file.is_open()){
      std::cout << " file not opened" << std::endl;
    }
    else
      {
	area_globalID.clear();
	area_cellArea.clear();
	while (!area_file.eof()) {
	  int globalID;
	  float cellarea;
	  area_file>>globalID>>cellarea;
	  if (area_file.eof()) break;
	  area_globalID.push_back(globalID);
	  area_cellArea.push_back(cellarea);	    
      }
	area_file.close();
      }
    std::cout<< "entries in the input text file "<<"\t"<<area_globalID.size()<<std::endl;
    //  Reading input calib files  //
    std::ifstream calib_file;
    calib_file.open(input_calib,ios::in);
    std::vector<int> calib_layer;
    std::vector<int> calib_globalID;
    std::vector<int> calib_RocChip;
    std::vector<int> calib_HalfRocChip;
    std::vector<int> calib_HalfRocChannel;
    std::vector<float> calib_gaussMean;
    std::vector<float> calib_gaussSigma;
    std::vector<float> calib_HistMean;
    int globalID_counter;
    cout<< "extracting pedestal fro module - "<<atoi(layer_)<<"\t"<<layer_<<endl;
    if(!calib_file.is_open()){
      std::cout << " file not opened" << std::endl;
    }
    else
      {
	calib_layer.clear();
	calib_globalID.clear();
	calib_RocChip.clear();
	calib_HalfRocChip.clear();
	calib_HalfRocChannel.clear();
	calib_gaussMean.clear();
	calib_gaussSigma.clear();
	calib_HistMean.clear();
	while (!calib_file.eof()) {
	  TString module_ ;
	  int module,layer,global_ID,chip_,halfchip,halfch,channel_,channel_map,channeltype,en_chan,entry;
	  float adc,adc_err,adc_sigma,err_sigma,chi2,grass, grass_err, mpv;
	  calib_file>>module_>>adc>>adc_sigma>>adc_err>>err_sigma>>mpv>>grass_err; //global_ID>>chip_>>halfchip>>halfch>>entry>>adc>>adc_sigma>>adc_err>>err_sigma>>mpv;
	  if (calib_file.eof()) break;
	  //cout<<module_<<"\t"<<atoi(layer_)<<"\t"<<layer<<endl;
	  if(module_.Contains("0x54")) layer = 0;
	  else  layer =1;
	  if(layer!=atoi(layer_)) continue;
	  //	  cout<<module_<<"\t"<<atoi(layer_)<<"\t"<<layer<<endl;

	  calib_layer.push_back(layer);
	  calib_globalID.push_back(globalID_counter);
	  // calib_RocChip.push_back(chip_);
	  // calib_HalfRocChip.push_back(halfchip);
	  // calib_HalfRocChannel.push_back(halfch);
	  
	  globalID_counter++;
	  
	  calib_gaussMean.push_back(adc);
	  calib_gaussSigma.push_back(adc_sigma);
	  //calib_HistMean.push_back(mpv);
	}
	calib_file.close();
      }
    std::cout<< "entries in the input text file "<<"\t"<<calib_gaussSigma.size()<<std::endl;
    
    // --------------------///
    // Reading pedestal input text files //
    // -------------------- ///
    std::ifstream in_file;
    in_file.open(inputfile1,ios::in);
    std::vector<int> map_layer;
    std::vector<int> map_globalID;
    std::vector<int> map_RocChip;
    std::vector<int> map_HalfRocChip;
    std::vector<int> map_HalfRocChannel;
    std::vector<float> map_gaussMean;
    std::vector<float> map_gaussSigma;
    std::vector<float> map_HistMean;

    std::vector<float> map_chi2;
    cout<< "extracting pedestal fro module - "<<atoi(layer_)<<"\t"<<layer_<<endl;
    if(!in_file.is_open())
      {
	std::cout << " file not opened" << std::endl;
      }
    else
      {
	 map_layer.clear();
	 map_globalID.clear();
	 map_RocChip.clear();
	 map_HalfRocChip.clear();
	 map_HalfRocChannel.clear();
	 map_gaussMean.clear();
	 map_gaussSigma.clear();
	 map_HistMean.clear();
	 while (!in_file.eof()) {
	   int layer,global_ID,chip_,halfchip,halfch,channel_,channel_map,channeltype,en_chan,entry;
	   float adc,adc_err,adc_sigma,err_sigma,chi2,grass, grass_err, mpv;
	   in_file>>layer>>global_ID>>chip_>>halfchip>>halfch>>entry>>adc>>adc_sigma>>adc_err>>err_sigma>>mpv>>chi2;
	   if (in_file.eof()) break;
	   if(layer!=atoi(layer_)) continue;
	   map_layer.push_back(layer);
	   map_globalID.push_back(global_ID);
	   map_RocChip.push_back(chip_);
	   map_HalfRocChip.push_back(halfchip);
	   map_HalfRocChannel.push_back(halfch);
	   map_gaussMean.push_back(adc);
	   map_gaussSigma.push_back(adc_sigma);
	   map_HistMean.push_back(mpv);
	   map_chi2.push_back(chi2);
	 }
	 in_file.close();
      }
    std::cout<< "entries in the input text file "<<"\t"<<map_HistMean.size()<<std::endl;

    int scheme = 0;
    TString title;
    //--------------------------------------------------
    // Test profile
    //--------------------------------------------------
    TProfile *profile = new TProfile("profile", "profile", 234, 0, 234, 0, 1024);
    switch(scheme) {
        default:
	  for(int i=0; i<234; ++i) {
	      float value = (float)map_HistMean[i] ;
	      if(ped_method==0)
		value = (float)map_HistMean[i];
	      else if(ped_method==1)
		value =(float)map_gaussMean[i];
	      else if(ped_method==2)
		value =(float)map_gaussSigma[i];
	      else if(ped_method==3)
                value = 100*(map_HistMean[i]-calib_gaussMean[i])/calib_gaussMean[i];
              else if(ped_method==4)
                value = 100*(map_gaussMean[i]-calib_gaussMean[i])/calib_gaussMean[i];
              else if(ped_method==5)
                value = (map_gaussSigma[i]-calib_gaussSigma[i])/calib_gaussSigma[i];
	      else if (ped_method==6)
		value = map_chi2[i];
	      else if(ped_method==7)
		value = map_gaussSigma[i]*126.46/area_cellArea[i]; // noise*d/A
	      //	      cout<< "value "<<value<<"\t"<<map_gaussSigma[i]<<"\t"<<area_cellArea[i]<<endl;
	      if(value==0)
		value = value+1e-6;
                if(i==0) profile->Fill(i, value+1e-6);
                else profile->Fill(i, value);		
            }
            break;

        case 1: // scheme: expected injected channels
            title = "Manual specification";
            for(int i=0; i<234; ++i) {
	      float value = map_HistMean[i]/calib_gaussMean[i];
	      if(ped_method==0)
                value = (float)map_HistMean[i];
              else if(ped_method==1)
                value =(float)map_gaussMean[i];
              else if(ped_method==2)
                value =(float)map_gaussSigma[i];
              else if(ped_method==3)
                value = 100*(map_HistMean[i]-calib_gaussMean[i])/calib_gaussMean[i];
              else if(ped_method==4)
                value = 100*(map_gaussMean[i]-calib_gaussMean[i])/calib_gaussMean[i];
              else if(ped_method==5)
                value = (map_gaussSigma[i]-calib_gaussSigma[i])/calib_gaussSigma[i];
              else if (ped_method==6)
                value = map_chi2[i];
              else if(ped_method==7)
                value = map_gaussSigma[i]/(0.1*area_cellArea[i]);

                if(i==0) profile->Fill(i, value+1e-6);
                else if(i==20) profile->Fill(i, value);
                else if(i==40) profile->Fill(i, value);
                else if(i==60) profile->Fill(i, value);
                else if(i==78) profile->Fill(i, value);
                else if(i==98) profile->Fill(i, value);
                else if(i==118) profile->Fill(i, value);
                else if(i==138) profile->Fill(i, value);
                else if(i==156) profile->Fill(i, value);
                else if(i==176) profile->Fill(i, value);
                else if(i==196) profile->Fill(i, value);
                else if(i==216) profile->Fill(i, value);
                else profile->Fill(i, -300.);
            }
            break;

        case 2: // scheme: results displayed on DQM GUI
            title = "DQM GUI (with readout sequence)";
            for(int i=0; i<234; ++i) {
	      float value = (float)map_HistMean[i];
	      if(ped_method==0)
		value =(float)map_HistMean[i];
	      else if(ped_method==1)
		value =(float)map_gaussMean[i];
	      else if(ped_method==2)
		value =(float)map_gaussSigma[i];

                if(i==21) profile->Fill(i, value);
                else if(i==42) profile->Fill(i, value);
                else if(i==64) profile->Fill(i, value);
                else if(i==77) profile->Fill(i, value);
                else if(i==99) profile->Fill(i, value);
                else if(i==120) profile->Fill(i, value);
                else if(i==142) profile->Fill(i, value);
                else if(i==155) profile->Fill(i, value);
                else if(i==177) profile->Fill(i, value);
                else if(i==198) profile->Fill(i, value);
                else if(i==220) profile->Fill(i, value);
                else profile->Fill(i, -300.);
            }
            break;

    }

    // for(int ibin=0; ibin<=235; ++ibin) {
    //     double value = profile->GetBinContent(ibin);
    //     printf("ibin = %d, value = %f\n", ibin, value);
    // }

    //--------------------------------------------------
    // Hexagonal plots
    //--------------------------------------------------
    profile->Draw();
    c1->SaveAs("test.root");
    TFile *f = TFile::Open(inputfile,"R");
    //    TString layer__ = layer_;
    title = "Module:"+layer_+", Run:"+ runNumber;//"LD wafer with global channel id (readout sequence)";
    TH2Poly *p = new TH2Poly("hexagonal histograms", title, -1*range, range, -1*range-2, range-2);
    p->SetStats(0);
    p->GetXaxis()->SetTitle("x (cm)");
    p->GetYaxis()->SetTitle("y (cm)");
    p->GetYaxis()->SetTitleOffset(1.1);
    
    title = "Module:"+layer_+" Run-"+ runNumber;//"LD wafer with HGCROC pin/chan";
    TH2Poly *p_pin = new TH2Poly("p_pin", title, -1*range, range, -1*range-2, range-2);
    p_pin->SetStats(0);
    p_pin->GetXaxis()->SetTitle("x (cm)");
    p_pin->GetYaxis()->SetTitle("y (cm)");
    p_pin->GetYaxis()->SetTitleOffset(1.1);

    title = "Module:"+layer_+" Run-"+ runNumber;//"LD wafer with Si cell pad Id";
    TH2Poly *p_sicell = new TH2Poly("p_sicell", title, -1*range, range, -1*range-2, range-2);
    p_sicell->SetStats(0);
    p_sicell->GetXaxis()->SetTitle("x (cm)");
    p_sicell->GetYaxis()->SetTitle("y (cm)");
    p_sicell->GetYaxis()->SetTitleOffset(1.1);

    int counter = 0;
    TGraph *gr;
    TKey *key;
    TIter nextkey(gDirectory->GetListOfKeys());
    while ((key = (TKey*)nextkey())) {
        TObject *obj = key->ReadObj();
        if(obj->InheritsFrom("TGraph")) {
            gr = (TGraph*) obj;
            p->AddBin(gr);
            p_pin->AddBin(gr);
            p_sicell->AddBin(gr);
            counter+=1;
        }
    }

    TRandom r;
    p->ChangePartition(100, 100);
    std::cout<<"counter "<<counter<<std::endl;
    for(int i=0; i<counter; ++i) {
        float value = profile->GetBinContent(i+1);
	if(value==0)
	  value = value+1e-06;
	//	if(i==146 || i==147)
	cout<< "value "<<value<<"\t"<<map_gaussSigma[i]<<"\t"<<area_cellArea[i]<<"\t"<<i<<endl;
        p->SetBinContent(i+1, value);

        // if(i==0 || i==78 || i==156)
        //     p_pin->SetBinContent(i+1, 0.000001);
        // else
	  p_pin->SetBinContent(i+1,value);// map_HGCROC_pin[i]);

        p_sicell->SetBinContent(i+1, value);//map_SiCell_padId[i]);
    }
float xmax=200;
    if(ped_method==0)
      {p->GetZaxis()->SetTitle("<ADC> counts");
	p_pin->GetZaxis()->SetTitle("<ADC> counts");
	p_sicell->GetZaxis()->SetTitle("<ADC> counts");
      }
    else if (ped_method==1)
      {p->GetZaxis()->SetTitle("#mu(ADC) -gauss");
        p_pin->GetZaxis()->SetTitle("#mu(ADC) -gauss");
        p_sicell->GetZaxis()->SetTitle("#mu (ADC) -gauss");
      }
    else if(ped_method==2)
      {
	p->GetZaxis()->SetTitle("#sigma(ADC) (gauss)");
        p_pin->GetZaxis()->SetTitle("#sigma(ADC) (gauss)");
        p_sicell->GetZaxis()->SetTitle("#sigma(ADC) (gauss)");
      xmax =5.0;
      }

    else if(ped_method==3)
      {p->GetZaxis()->SetTitle("100*(<ADC> - <Calib ADC >)/<Calib ADC>");
        p_pin->GetZaxis()->SetTitle("100*(<ADC> - <Calib ADC >)/<Calib ADC>");
        p_sicell->GetZaxis()->SetTitle("100*(<ADC> - <Calib ADC >)/<Calib ADC>");
	xmax=1.5;
      }
    else if (ped_method==4)
      {p->GetZaxis()->SetTitle("100*(#mu- <Calib ADC >)/<Calib ADC>");
	p_pin->GetZaxis()->SetTitle("100*(#mu- <Calib ADC >)/<Calib ADC>");
        p_sicell->GetZaxis()->SetTitle("100*(#mu- <Calib ADC >)/<Calib ADC>");
	xmax =1.5;
      }
    else if(ped_method==5)
      {
        p->GetZaxis()->SetTitle("(#sigma- <Calib sigma >)/<Calib sigma>");
        p_pin->GetZaxis()->SetTitle("(#sigma- <Calib sigma >)/<Calib sigma>");
        p_sicell->GetZaxis()->SetTitle("(#sigma- <Calib sigma >)/<Calib sigma>");
	xmax =1.1;
      }

    else if(ped_method==6)
      {
        p->GetZaxis()->SetTitle("#chi^{2}/ndf");
        p_pin->GetZaxis()->SetTitle("#chi^{2}/ndf");
        p_sicell->GetZaxis()->SetTitle("#chi^{2}/ndf");
	xmax =500;
      }

    else if(ped_method==7)
      {
        p->GetZaxis()->SetTitle("Noise (Normalized to full cell area)");
        p_pin->GetZaxis()->SetTitle("Noise (Normalized to full cell area)");
        p_sicell->GetZaxis()->SetTitle("Noise (Normalized to full cell area)");
        xmax =5;
      }

    p->SetMarkerSize(0.8);
    p_pin->SetMarkerSize(0.8);
    p_sicell->SetMarkerSize(0.8);

    p->SetMaximum(xmax);
    p->GetZaxis()->SetTitleOffset(1.3);
    p_pin->GetZaxis()->SetTitleOffset(1.3);
    p_sicell->GetZaxis()->SetTitleOffset(1.3);

    p->Draw("colz;text");
    beautify_plot();
    //    c1->SaveAs(outputfile);
    TString  eosdir = "/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Aug2023_TBHGCAL_CERN/PedestalStudies/Hexplots/";
    c1->SaveAs(eosdir+png_String+"_globalChannelId_readoutSequence.png");
    c1->SaveAs(eosdir+png_String+"_globalChannelId_readoutSequence.pdf");
    c1->SaveAs(png_String+"_globalChannelId_readoutSequence.png");
    c1->SaveAs(png_String+"_globalChannelId_readoutSequence.pdf");

    p_pin->SetMarkerSize(0.7);
    p_pin->SetMaximum(xmax);

    p_pin->Draw("colz;text");
    beautify_plot();

    c1->SaveAs(png_String+"_HGCROC_pin_chan.png");
    c1->SaveAs(png_String+"_HGCROC_pin_chan.pdf");
    
    p_sicell->SetMaximum(xmax);
    p_sicell->SetMarkerSize(0.7);
    p_sicell->Draw("colz;text");
    beautify_plot();
    c1->SaveAs(png_String+"_SiCell_padId.png");
    c1->SaveAs(png_String+"_SiCell_padId.pdf");

    //-----------------------------------------------------------------
    // Reminder: counter = nCells - 9
    //-----------------------------------------------------------------
    // printf("[INFO] nCells  = %d\n", p->GetNcells());
    // printf("[INFO] counter = %d\n", counter);

}

void beautify_plot(bool drawLine = true, bool drawText = true) {
    //-----------------------------------------------------------------
    // cosmetics
    //-----------------------------------------------------------------
    if(drawLine) {
        // load N_boundary_points, x1, x2, x3, y1, y2, y3 from auxiliary_boundary_lines.h

        TLine line;
        line.SetLineStyle(1);
        line.SetLineColor(2);
        line.SetLineWidth(2);

        for(int i=0; i<aux::N_boundary_points-1; ++i) {
            line.DrawLine(aux::x1[i], aux::y1[i], aux::x1[i+1], aux::y1[i+1]);
            line.DrawLine(aux::x2[i], aux::y2[i], aux::x2[i+1], aux::y2[i+1]);
            line.DrawLine(aux::x3[i], aux::y3[i], aux::x3[i+1], aux::y3[i+1]);
            line.DrawLine(aux::x4[i], aux::y4[i], aux::x4[i+1], aux::y4[i+1]);
            line.DrawLine(aux::x5[i], aux::y5[i], aux::x5[i+1], aux::y5[i+1]);
            line.DrawLine(aux::x6[i], aux::y6[i], aux::x6[i+1], aux::y6[i+1]);
        }
    }

    if(drawText) {
        TText text;
        text.SetTextAlign(22);
        text.SetTextFont(43);
        text.SetTextSize(12);

        double theta1 = -TMath::Pi()/3.;
        double theta2 = TMath::Pi()/3.;
        double theta3 = TMath::Pi();
        std::vector<double> theta_angle_text = {60, 60, -60, -60, 0, 0};
        std::vector<double> theta_coordinate_text = {theta1, theta1, theta2, theta2, theta3, theta3};
        //std::vector<double> x_coordinate_text = {-6.25, 6.25, -6.25, 6.25, -6.25, 6.25};
        std::vector<double> x_coordinate_text = {-6.25, 6.25, -5, 7.5, -6.25, 6.25};
        std::vector<double> y_coordinate_text = {23.5, 23.5, 26, 26, 26.4, 26.4};
        std::vector<TString> v_texts = {"chip-0, half-1", "chip-0, half-0",
                                        "chip-1, half-1", "chip-1, half-0",
                                        "chip-2, half-1", "chip-2, half-0"};

        double arbUnit_to_cm = 17./24.;

        // evaluate (r, phi) and apply rotation
        for(int i=0; i<6; ++i) {
            text.SetTextAngle(theta_angle_text[i]);
            double theta = theta_coordinate_text[i];
            double cos_theta = TMath::Cos(theta);
            double sin_theta = TMath::Sin(theta);

            double x = x_coordinate_text[i];
            double y = y_coordinate_text[i];
            double r = sqrt(pow(x,2)+pow(y,2));
            double cos_phi = x/r;
            double sin_phi = y/r;
            x = r*(cos_phi*cos_theta + sin_phi*sin_theta)*arbUnit_to_cm;
            y = r*(sin_phi*cos_theta - cos_phi*sin_theta)*arbUnit_to_cm;
            text.DrawText(x, y, v_texts[i]);
        }
    }
}
