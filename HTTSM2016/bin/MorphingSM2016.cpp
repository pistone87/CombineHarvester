#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include "boost/algorithm/string/predicate.hpp"
#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/HTTSM2016/interface/HttSystematics_SMRun2.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TH2.h"
#include "TF1.h"

using namespace std;
using boost::starts_with;
namespace po = boost::program_options;

template <typename T>
void To1Bin(T* proc)
{
    std::unique_ptr<TH1> originalHist = proc->ClonedScaledShape();
    TH1F *hist = new TH1F("hist","hist",1,0,1);
    double err = 0;
    double rate =
    originalHist->IntegralAndError(0, originalHist->GetNbinsX() + 1, err);
    hist->SetDirectory(0);
    hist->SetBinContent(1, rate);
    hist->SetBinError(1, err);
    proc->set_shape(*hist, true);  // True means adjust the process rate to the integral of the hist
}

bool BinIsControlRegion(ch::Object const* obj)
{
    return (boost::regex_search(obj->bin(),boost::regex{"_cr$"}) || (obj->channel() == std::string("mm")));
}

// Useful to have the inverse sometimes too
bool BinIsNotControlRegion(ch::Object const* obj)
{
    return !BinIsControlRegion(obj);
}



int main(int argc, char** argv) {

    string SM125= "";
    string mass = "mA";
    string output_folder = "sm_run2";
    string only_init = "";
    string postfix="";
	string scale_sig_procs="";

    string input_folder_em="USCMS/";
    string input_folder_et="USCMS/";
    string input_folder_mt="USCMS/";
    string input_folder_tt="USCMS/";
    string input_folder_mm="USCMS/";
    string input_folder_ttbar="USCMS/";

    int control_region = 0;
    bool auto_rebin = false;
    bool manual_rebin = false;
    bool real_data = false;
    bool check_neg_bins = false;
    bool poisson_bbb = false;
    bool do_w_weighting = false;
    bool mm_fit = false;
    bool ttbar_fit = false;
    bool do_jetfakes = false;

    po::variables_map vm;
    po::options_description config("configuration");
    config.add_options()

    ("SM125,h", po::value<string>(&SM125)->default_value(SM125))
    ("mass,m", po::value<string>(&mass)->default_value(mass))
    ("output_folder", po::value<string>(&output_folder)->default_value("sm_run2"))
    ("only_init", po::value<string>(&only_init)->default_value(""))
    ("postfix", po::value<string>(&postfix)->default_value(""))
	("scale_sig_procs", po::value<string>(&scale_sig_procs)->default_value(""))

    ("input_folder_em", po::value<string>(&input_folder_em)->default_value("USCMS"))
    ("input_folder_et", po::value<string>(&input_folder_et)->default_value("USCMS"))
    ("input_folder_mt", po::value<string>(&input_folder_mt)->default_value("USCMS"))
    ("input_folder_tt", po::value<string>(&input_folder_tt)->default_value("USCMS"))
    ("input_folder_mm", po::value<string>(&input_folder_mm)->default_value("USCMS"))
    ("input_folder_ttbar", po::value<string>(&input_folder_ttbar)->default_value("USCMS"))

    ("control_region", po::value<int>(&control_region)->default_value(0))
    ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(false))
    ("manual_rebin", po::value<bool>(&manual_rebin)->default_value(false))
    ("real_data", po::value<bool>(&real_data)->default_value(false))
    ("check_neg_bins", po::value<bool>(&check_neg_bins)->default_value(false))
    ("poisson_bbb", po::value<bool>(&poisson_bbb)->default_value(false))
    ("w_weighting", po::value<bool>(&do_w_weighting)->default_value(false))
    ("mm_fit", po::value<bool>(&mm_fit)->default_value(true))
    ("ttbar_fit", po::value<bool>(&ttbar_fit)->default_value(true))
    ("jetfakes", po::value<bool>(&do_jetfakes)->default_value(false));

    po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
    po::notify(vm);



    typedef vector<string> VString;
    typedef vector<pair<int, string>> Categories;
    //! [part1]
    // First define the location of the "auxiliaries" directory where we can
    // source the input files containing the datacard shapes
    //    string aux_shapes = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/CombineTools/bin/AllROOT_20fb/";
    std::map<string, string> input_dir;
    input_dir["em"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSM2016/shapes/"+input_folder_em+"/";
    input_dir["mt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSM2016/shapes/"+input_folder_mt+"/";
    input_dir["et"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSM2016/shapes/"+input_folder_et+"/";
    input_dir["tt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSM2016/shapes/"+input_folder_tt+"/";
    input_dir["mm"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSM2016/shapes/"+input_folder_mm+"/";
    input_dir["ttbar"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSM2016/shapes/"+input_folder_ttbar+"/";


    VString chns = {"mt","et","tt","em"};
    if (mm_fit) chns.push_back("mm");
    if (ttbar_fit) chns.push_back("ttbar");


    map<string, VString> bkg_procs;
    if (do_jetfakes){
      bkg_procs["et"] = {"ZTT", "ZL", "TTT", "VVT", "EWKZ", "jetFakes"};
      bkg_procs["mt"] = {"ZTT", "ZL", "TTT", "VVT", "EWKZ", "jetFakes"};
      bkg_procs["tt"] = {"ZTT", "ZL", "TTT", "VVT", "EWKZ", "jetFakes", "W_rest", "ZJ_rest", "TTJ_rest","VVJ_rest"};
    }else{
      bkg_procs["et"] = {"ZTT", "QCD", "ZL", "ZJ","TTT","TTJ", "VV", "EWKZ"};
      bkg_procs["mt"] = {"ZTT", "QCD", "ZL", "ZJ","TTT","TTJ", "VV", "EWKZ"};
      bkg_procs["tt"] = {"ZTT", "W", "QCD", "ZL", "ZJ","TTT","TTJ", "VVT", "VVJ", "EWKZ"};
    }
    bkg_procs["em"] = {"ZTT", "W", "QCD", "ZL", "TT", "VV", "EWKZ", "ggH_hww125", "qqH_hww125"};
    bkg_procs["mm"] = {"W", "ZL", "TT", "VV"};
    bkg_procs["ttbar"] = {"ZTT", "W", "QCD", "ZL", "TT", "VV", "EWKZ"};


    ch::CombineHarvester cb;


    map<string,Categories> cats;
	// cats for et channel
    cats["et"] = {
        {1, "et_0jet"},
        {2, "et_boosted"},
        {3, "et_vbf"},
    };
	// cats for mt channel
    cats["mt"] = {
        {1, "mt_0jet"},
        {2, "mt_boosted"},
        {3, "mt_vbf"},
    };
	// cats for em channel
    cats["em"] = {
        {1, "em_0jet"},
        {2, "em_boosted"},
        {3, "em_vbf"}
	};
	// cats for tt channel
    cats["tt"] = {
        {1, "tt_0jet"},
        {2, "tt_boosted"},
        {3, "tt_vbf"}
	};
	// cats for mm channel
    cats["mm"] = {
        {1, "mm_0jet"},
        {2, "mm_1jet"},
        {3, "mm_vbf"}
    };
	// cats for tt channel
    cats["ttbar"] = {
        {1, "ttbar_all"}
    };


    if (control_region > 0){
        // for each channel use the categories >= 10 for the control regions
        // the control regions are ordered in triples (10,11,12),(13,14,15)...
        for (auto chn : chns){
            // for em or mm do nothing

			// et and mt channels
            if (ch::contains({"mt", "et"}, chn)) {
                Categories queue;
                int binid = 10;
                queue.push_back(make_pair(binid,chn+"_wjets_0jet_cr"));
                queue.push_back(make_pair(binid+1,chn+"_wjets_boosted_cr"));
                //queue.push_back(make_pair(binid+2,chn+"_wjets_vbf_cr"));
                queue.push_back(make_pair(binid+3,chn+"_antiiso_0jet_cr"));
                queue.push_back(make_pair(binid+4,chn+"_antiiso_boosted_cr"));
                //queue.push_back(make_pair(binid+5,chn+"_antiiso_vbf_cr"));

                cats[chn].insert(cats[chn].end(),queue.begin(),queue.end());
            } // et or mt channel
            // tt channel
            if (ch::contains({"tt"}, chn)) {
                Categories queue;
                int binid = 10;
                queue.push_back(make_pair(binid,chn+"_0jet_qcd_cr"));
                queue.push_back(make_pair(binid+1,chn+"_boosted_qcd_cr"));
                queue.push_back(make_pair(binid+2,chn+"_vbf_qcd_cr"));

                cats[chn].insert(cats[chn].end(),queue.begin(),queue.end());
            } // tt channel
        } // end loop over channels
    } // end control_region > 0




    // Or equivalently, specify the mass points explicitly
    //vector<string> sig_procs = {"ggH_htt","qqH_htt","WH_htt","ZH_htt"};
	// the processes names are the same given for the systematics
	// "ALT" is needed for the hypothesis test
	vector<string> sig_procs = {"smHcpeven", "susyHcpodd_ALT"};
    vector<string> masses = {"125"};

    using ch::syst::bin_id;

    //! [part2]
    for (auto chn : chns) {
        cb.AddObservations({"*"}, {"htt"}, {"13TeV"}, {chn}, cats[chn]);
        cb.AddProcesses(   {"*"}, {"htt"}, {"13TeV"}, {chn}, bkg_procs[chn], cats[chn], false);
        cb.AddProcesses(masses,   {"htt"}, {"13TeV"}, {chn}, sig_procs, cats[chn], true);
        //Needed to add ewkz and W as these are not not available/Negative in qcd cR
    }

    //Add EWKZ and W manually !!!!!!
    if (control_region > 0){
      cb.AddProcesses(   {"*"}, {"htt"}, {"13TeV"}, {"et"}, {"W"}, {{1, "et_0jet"},{2, "et_boosted"},{3, "et_vbf"},
                                    {10, "et_wjets_0jet_cr"},{11, "et_wjets_boosted_cr"},
                                    {13, "et_antiiso_0jet_cr"},{14, "et_antiiso_boosted_cr"}}, false);
      cb.AddProcesses(   {"*"}, {"htt"}, {"13TeV"}, {"mt"}, {"W"}, {{1, "mt_0jet"},{2, "mt_boosted"},{3, "mt_vbf"},
                                    {10, "mt_wjets_0jet_cr"},{11, "mt_wjets_boosted_cr"},
                                    {13, "mt_antiiso_0jet_cr"},{14, "mt_antiiso_boosted_cr"}}, false);
    }else if (!do_jetfakes){
      cb.AddProcesses(   {"*"}, {"htt"}, {"13TeV"}, {"et"}, {"W"}, {{1, "et_0jet"},{2, "et_boosted"},{3, "et_vbf"}}, false);
      cb.AddProcesses(   {"*"}, {"htt"}, {"13TeV"}, {"mt"}, {"W"}, {{1, "mt_0jet"},{2, "mt_boosted"},{3, "mt_vbf"}}, false);
    }



    //! [part4]
    if ((control_region > 0) || mm_fit){
        // Since we now account for QCD in the high mT region we only
        // need to filter signal processes
        cb.FilterAll([](ch::Object const* obj) {
            return (BinIsControlRegion(obj) && obj->signal());
        });
    }


    ch::AddSMRun2Systematics(cb, control_region, mm_fit, ttbar_fit);

    if (! only_init.empty()) {
        std::cout << "Write datacards (without shapes) to directory \"" << only_init << "\" and quit." << std::endl;
        ch::CardWriter tmpWriter("$TAG/$ANALYSIS_$ERA_$CHANNEL_$BINID_$MASS.txt", "$TAG/dummy.root");
        tmpWriter.WriteCards(only_init, cb);
        return 0;
    }

    //! [part7]
    //for (string chn:chns){
	for (string chn : cb.channel_set()){
		string channel = chn;
        cb.cp().channel({chn}).backgrounds().ExtractShapes(
                                                           input_dir[chn] + "htt_"+chn+".inputs-sm-13TeV"+postfix+".root",
                                                           "$BIN/$PROCESS",
                                                           "$BIN/$PROCESS_$SYSTEMATIC");
        cb.cp().channel({chn}).process(sig_procs).ExtractShapes(
                                                                input_dir[chn] + "htt_"+chn+".inputs-sm-13TeV"+postfix+".root",
                                                                "$BIN/$PROCESS$MASS",
                                                                "$BIN/$PROCESS$MASS_$SYSTEMATIC");
    }




    // Now delete processes with 0 yield
    cb.FilterProcs([&](ch::Process *p) {
        bool null_yield = !(p->rate() > 0. || BinIsControlRegion(p));
        if (null_yield){
            std::cout << "[Null yield] Removing process with null yield: \n ";
            std::cout << ch::Process::PrintHeader << *p << "\n";
            cb.FilterSysts([&](ch::Systematic *s){
                bool remove_syst = (MatchingProcess(*p,*s));
                return remove_syst;
            });
        }
        return null_yield;
    });



    // And convert any shapes in the CRs to lnN:
    // Convert all shapes to lnN at this stage
    cb.cp().FilterSysts(BinIsNotControlRegion).syst_type({"shape"}).ForEachSyst([](ch::Systematic *sys) {
        sys->set_type("lnN");
		if(sys->value_d() <0.001) {sys->set_value_d(0.001);};
    });




    //! [part8]
    auto bbb = ch::BinByBinFactory()
    .SetAddThreshold(0.05)
    .SetMergeThreshold(0.8)
    .SetFixNorm(false);
    bbb.MergeBinErrors(cb.cp().backgrounds());
    bbb.AddBinByBin(cb.cp().backgrounds(), cb);



    // And now do bbb for the control region with a slightly different config
    auto bbb_ctl = ch::BinByBinFactory()
    .SetPattern("CMS_$ANALYSIS_$BIN_$ERA_$PROCESS_bin_$#")
    .SetAddThreshold(0.)
    .SetMergeThreshold(0.8)
    .SetFixNorm(false)  // contrary to signal region, bbb *should* change yield here
    .SetVerbosity(1);
    // Will merge but only for non W and QCD processes, to be on the safe side
    bbb_ctl.MergeBinErrors(cb.cp().process({"QCD", "W"}, false).FilterProcs(BinIsNotControlRegion));
    bbb_ctl.AddBinByBin(cb.cp().process({"QCD", "W"}, false).FilterProcs(BinIsNotControlRegion), cb);
    cout << " done\n";




    // This function modifies every entry to have a standardised bin name of
    // the form: {analysis}_{channel}_{bin_id}_{era}
    // which is commonly used in the htt analyses
    ch::SetStandardBinNames(cb);


    // Write out datacards. Naming convention important for rest of workflow. We
    // make one directory per chn-cat, one per chn and cmb. In this code we only
    // store the individual datacards for each directory to be combined later, but
    // note that it's also possible to write out the full combined card with CH
    string output_prefix = "output/";
    if(output_folder.compare(0,1,"/") == 0) output_prefix="";
    ch::CardWriter writer(output_prefix + output_folder + "/$TAG/$MASS/$BIN.txt",
                          output_prefix + output_folder + "/$TAG/common/htt_input.root");


    writer.WriteCards("cmb", cb);
    for (auto chn : chns) {
        if(chn == std::string("mm"))
        {
            continue;
        }
        // per-channel
        writer.WriteCards(chn, cb.cp().channel({chn, "mm"}));
        // And per-channel-category
        writer.WriteCards("htt_"+chn+"_1_13TeV", cb.cp().channel({chn}).bin_id({1}));
        writer.WriteCards("htt_"+chn+"_2_13TeV", cb.cp().channel({chn}).bin_id({2}));
        writer.WriteCards("htt_"+chn+"_3_13TeV", cb.cp().channel({chn}).bin_id({3}));


        for (auto mmm : {"125"}){

            if (mm_fit){
                cb.cp().channel({"mm"}).bin_id({1}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_mm_1_13TeV.txt", output_prefix + output_folder +"/"+chn+"/common/htt_input_mm1.root");
                cb.cp().channel({"mm"}).bin_id({2}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_mm_2_13TeV.txt", output_prefix + output_folder +"/"+chn+"/common/htt_input_mm2.root");
                cb.cp().channel({"mm"}).bin_id({3}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_mm_3_13TeV.txt", output_prefix + output_folder +"/"+chn+"/common/htt_input_mm3.root");
            }

            if (ttbar_fit){
                cb.cp().channel({"ttbar"}).bin_id({1}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_ttbar_1_13TeV.txt", output_prefix + output_folder +"/"+chn+"/common/htt_input_ttbar1.root");
            }


            if (control_region > 0){
				// et or mt channels
                if (ch::contains({"et", "mt"}, chn)) {

                    cb.cp().channel({chn}).bin_id({10}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_"+chn+"_10_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"10.root");
                    cb.cp().channel({chn}).bin_id({11}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_"+chn+"_11_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"11.root");
                    //cb.cp().channel({chn}).bin_id({12}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_"+chn+"_12_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"12.root");
                    cb.cp().channel({chn}).bin_id({13}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_"+chn+"_13_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"13.root");
                    cb.cp().channel({chn}).bin_id({14}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_"+chn+"_14_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"14.root");
                    //cb.cp().channel({chn}).bin_id({15}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_"+chn+"_15_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"15.root");


                    cb.cp().channel({chn}).bin_id({10}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/cmb/"+mmm+ "/htt_"+chn+"_10_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"10.root");
                    cb.cp().channel({chn}).bin_id({11}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/cmb/"+mmm+ "/htt_"+chn+"_11_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"11.root");
                    //cb.cp().channel({chn}).bin_id({12}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/cmb/"+mmm+ "/htt_"+chn+"_12_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"12.root");
                    cb.cp().channel({chn}).bin_id({13}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/cmb/"+mmm+ "/htt_"+chn+"_13_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"13.root");
                    cb.cp().channel({chn}).bin_id({14}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/cmb/"+mmm+ "/htt_"+chn+"_14_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"14.root");
                    //cb.cp().channel({chn}).bin_id({15}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/cmb/"+mmm+ "/htt_"+chn+"_15_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"15.root");

                } // end if et or mt
				// tt channel
                if (ch::contains({"tt"}, chn)) {

                    cb.cp().channel({chn}).bin_id({10}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_"+chn+"_10_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"10.root");
                    cb.cp().channel({chn}).bin_id({11}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_"+chn+"_11_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"11.root");
                    cb.cp().channel({chn}).bin_id({12}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/"+chn+"/"+mmm+ "/htt_"+chn+"_12_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"12.root");


                    cb.cp().channel({chn}).bin_id({10}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/cmb/"+mmm+ "/htt_"+chn+"_10_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"10.root");
                    cb.cp().channel({chn}).bin_id({11}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/cmb/"+mmm+ "/htt_"+chn+"_11_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"11.root");
                    cb.cp().channel({chn}).bin_id({12}).mass({"$MASS", "*"}).WriteDatacard(output_prefix + output_folder +"/cmb/"+mmm+ "/htt_"+chn+"_12_13TeV.txt", output_prefix + output_folder +"/"+chn+ "/common/htt_input"+chn+"12.root");
                } // end f tt channel
            } // end control_region
        }
    }


    cb.PrintAll();
    cout << " done\n";


}
