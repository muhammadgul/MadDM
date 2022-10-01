    #include <iostream>
    #include <fstream>
    #include <string>
    #include <map>
    #include "TH1.h"
    #ifdef __MAKECINT__
    #pragma link C++ class std::map<std::string,TH1*>+;
    #endif
    
    void same_canvas()
    { bool mass=0;
      TCanvas*c1 = new TCanvas("c1","E_{T}^{miss} for tan#{beta}",600,600);
     
      TFile* fIn1 = new TFile("lhe_files/for_fig14_left_13TeV/run_09/unweighted_events.root");
      TFile* fIn5 = new TFile("lhe_files/for_fig14_left_13_6TeV/run_04/unweighted_events.root");
      TFile* fIn10 = new TFile("lhe_files/for_fig14_left_27TeV/run_01/unweighted_events.root");
      TFile* fIn100 = new TFile("lhe_files/for_fig14_left_100TeV/run_01/unweighted_events.root");
           
      TH1F *H1 = (TH1F*)fIn1->Get("pt_xd_xd");
      TH1F *H2 = (TH1F*)fIn5->Get("pt_xd_xd");
      TH1F *H3 = (TH1F*)fIn10->Get("pt_xd_xd");
      TH1F *H4 = (TH1F*)fIn100->Get("pt_xd_xd");
      TH1F *H1clone =  (TH1F*)  H1->Clone();
      TH1F *H2clone =  (TH1F*)  H2->Clone();
      TH1F *H3clone =  (TH1F*)  H3->Clone();
      TH1F *H4clone =  (TH1F*)  H4->Clone();
  cout<<"these are bins: "<<H1->GetNbinsX()<<endl;
/*      for (int i=1; i<= H1->GetNbinsX(); i++){
        double ratio1 = H1->GetBinContent(i)/(H1->GetBinWidth(i)*H1->Integral());
        double ratio2 = H2->GetBinContent(i)/(H2->GetBinWidth(i)*H2->Integral());
        double ratio3 = H3->GetBinContent(i)/(H3->GetBinWidth(i)*H3->Integral());
        H1clone->SetBinContent(i,ratio1); 
        H2clone->SetBinContent(i,ratio2); 
        H3clone->SetBinContent(i,ratio3); 
     */
    /*    H1->Draw("hist");
        H1->GetYaxis()->SetRangeUser(-1000.0,5000.);
 
        H2->Draw("hist same"); 
        H3->Draw("hist same"); 
      */  
//        cout<<"integ: "<<H1clone->Integral()<<endl;
//  H1->Scale(1/H1->Integral());   
//  H2->Scale(1/H2->Integral());   
//  H3->Scale(1/H3->Integral());   
      H1clone->SetLineColor(kGreen);
      H1clone->SetLineWidth(3);
      H2clone->SetLineColor(kBlue);
      H2clone->SetLineWidth(3);
      H3clone->SetLineColor(kRed);
      H3clone->SetLineWidth(3);
      H4clone->SetLineColor(kMagenta);
      H4clone->SetLineWidth(3);
// tanb, xsec (1, .0028536 pb) (5, 0.0718 pb) (10, 0.2847 pb)
//     H2->SetMaximum(H2->GetMaximum()*520000);
//     H3clone->DrawNormalized("hist",(0.2847));
cout<<"integral: "<<H3clone->Integral()<<endl;
     H3clone->Scale(1/H3clone->Integral());
//     H1clone->DrawNormalized("hist same",(0.0028536));
     H1clone->Scale(1/H1clone->Integral());
//     H2clone->DrawNormalized("hist same",(0.0718));
     H2clone->Scale(1/H2clone->Integral());
     H4clone->Scale(1/H4clone->Integral());

      H4clone->Draw("hist");
      H3clone->Draw("hist same");
      H2clone->Draw("hist same");
      H1clone->Draw("hist same");
//     H1clone->GetYaxis()->SetRangeUser(0.000004,1.50);


//      H1->GetYaxis()->SetRange(0.,50.);

       //    H2->SetMarkerStyle(21);
//      H1->Draw("hist");
//      H3clone->SetMaximum(H3clone->GetMaximum()*5.5);
//      H2->Draw("hist same");
//      H3->Draw("hist same");
  //     TGaxis::SetMaxDigits(2);
//      H1clone->SetMarkerStyle(21);
      H4clone->SetTitle("");
      H4clone->GetXaxis()->SetTitle("#it{E}^{miss}_{#it{T}}");
      H4clone->GetYaxis()->SetTitle("Normalized");
      H4clone->GetYaxis()->SetTitleOffset(1.50);
      H4clone->SetStats(0);

    TLatex *tex1 = new TLatex(0.10,0.92,"mono-Higgs  gg fusion");
    tex1->SetNDC();
    tex1->SetTextAngle(0);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04);
    tex1->SetTextAlign(11);
    tex1->Draw();
    TLatex *tex2 = new TLatex(0.565,0.92,"M_{A,a} = {700, 200 GeV}");
    tex2->SetNDC();
    tex2->SetTextAngle(0);
    tex2->SetTextFont(42);
    tex2->SetTextAlign(11);
    tex2->SetTextSize(0.04);
    tex2->Draw();
    TLatex *tex3 = new TLatex(0.45,0.8,"tab#beta = 1");
    tex3->SetNDC();
    tex3->SetTextAngle(0);
    tex3->SetTextFont(42);
    tex3->SetTextAlign(11);
    tex3->SetTextSize(0.04);
    tex3->Draw();

 //     pdfHistclone[0]->GetYaxis()->SetRangeUser(0.92,1.08);
  //      TLegend *leg = new TLegend(0.15,0.75,0.55,0.85,NULL,"brNDC");//1)move to left and expand 2)move down and vertically seperate 3)move left and shrink 4) shrink vertically
        TLegend *leg = new TLegend(0.15,0.6,0.60,0.8,NULL,"brNDC");//1)move to left and expand 2)move down and vertically seperate 3)move left and shrink 4) shrink vertically
//        leg-> SetNColumns(3);
//        leg->SetFillStyle ( 0);
//        leg->SetFillColor ( 0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(10);
        leg->SetTextSize(0.04);
        leg->Draw();
        leg->AddEntry(H1clone,"#sqrt{s} = 13 TeV","LP");
        leg->AddEntry(H2clone,"#sqrt{s} = 13.6 TeV","LP");
        leg->AddEntry(H3clone,"#sqrt{s} = 27 TeV","LP");
        leg->AddEntry(H4clone,"#sqrt{s} = 100 TeV","LP");
        c1->Update();
//        c1->SetLogy();
//        c1->SetGrid();
//        c1->SetTicks();
        c1->SaveAs("results/tan_beta1_13_vs13_6_vs_27_vs_100TeV.pdf");
      
    }
