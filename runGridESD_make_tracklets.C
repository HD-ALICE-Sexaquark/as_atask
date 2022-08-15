class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisManager;
class AliPhysicsSelectionTask;
class AliCentralitySelectionTask;
class AliEmcalCorrectionTask;
class AliEmcalJetTask;
class AliAnalysisGrid;
class AliAnalysisAlien;

#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "Ali_make_tracklets_from_digits.h"

void runGridESD_make_tracklets() {
    // aliroot runGridESD_make_tracklets.C
    gInterpreter->GenerateDictionary("std::vector<std::vector<std::vector<std::vector<Double_t>>>>", "vector");

    TString fname = "TRD_tracklets";

    Int_t beamtime = 7;  // 0: p-Pb 2016, 1: Pb-Pb 2018, 2: p-p simulation, 6: p-Pb nuclei simulatio

    Bool_t isMC = kFALSE;   // Monte Carlo or "real" data
    Bool_t local = kFALSE;  // kTRUE for local analysis, kFALSE for grid analysis
    Int_t sub = 702;        // 702
    TString sRunPeriod;

    if (beamtime == 0)  // p-Pb 2016
    {
        sub = 702;
        sRunPeriod = "LHC16q";
    }
    if (beamtime == 1)  // Pb-Pb 2018
    {
        sub = 800;
        sRunPeriod = "LHC18r";
    }
    if (beamtime == 2)  // p-p simulation
    {
        sub = 900;
        isMC = kTRUE;
        sRunPeriod = "LHC18r";
    }
    if (beamtime == 3)  // p-Pb simulation
    {
        sub = 910;
        isMC = kTRUE;
        sRunPeriod = "LHC16q";
    }
    if (beamtime == 4)  // p-Pb 2016, one full run
    {
        sub = 7021;
        sRunPeriod = "LHC16q";
    }
    if (beamtime == 5)  // Krypton data
    {
        sub = 888;
        sRunPeriod = "LHC18r";
    }
    if (beamtime == 6)  // p-Pb nuclei simulation
    {
        sub = 123;
        isMC = kTRUE;
        sRunPeriod = "LHC22a3";
    }
    if (beamtime == 7)  // S-Particle
    {
        sub = 555;
        isMC = kTRUE;
        local = kTRUE;
        sRunPeriod = "LHC22a3";
    }

    /// set parameters for the analysis
    const Char_t *cDataType = "ESD";               // set analysis type; AOD or ESD
    const Char_t *cRunPeriod = sRunPeriod.Data();  // set run period, LHC18q
    const Char_t *cGridMode = "full";              // grid mode; test, full or terminate (for merging)
    Bool_t useJDL = kTRUE;
    const Char_t *cTaskName = "TRD_Make_Tracklets";  // name of the task

    /// since we will compile a class, tell root where to look for headers
#if !defined(__CINT__) || defined(__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
#endif

    /// get beam type (default: pp)
    // TString sRunPeriod(cRunPeriod);

    /// create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager(cTaskName);
    AliESDInputHandler *esdH = new AliESDInputHandler();
    esdH->SetFriendFileName("AliESDfriends.root");
    esdH->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdH);

    // AliAnalysisManager *mgr = new AliAnalysisManager(cTaskName);
    // AliAODInputHandler *aodH = new AliAODInputHandler();
    // mgr->SetInputEventHandler(aodH);

    /// centrality selection (multiplicity task)
    // TMacro multSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
    // AliMultSelectionTask* multSelectionTask = reinterpret_cast<AliMultSelectionTask*>(multSelection.Exec());

    // Centrality
    // gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    // AddTaskMultSelection();

    //               MULT SELECTION
    //#if !defined (__CINT__) || defined (__CLING__)
    // AliMultSelectionTask
    // *multSelTask=reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C()"));
    //#else
    // gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    // AddTaskMultSelection();
    //#endif

    // load the macro and add the task
    TMacro PIDadd(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
    AliAnalysisTaskPIDResponse *PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse *>(PIDadd.Exec());

    // gInterpreter->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    // gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    // AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(kFALSE,kTRUE,kTRUE,"1");

    // compile the class (locally) with debug symbols
    gInterpreter->LoadMacro("Ali_make_tracklets_from_digits.cxx+g");

    // load the addtask macro and create the task
    Ali_make_tracklets_from_digits *myTask = 0;
    myTask = reinterpret_cast<Ali_make_tracklets_from_digits *>(gInterpreter->ExecuteMacro(Form("AddTask_tracklets_aschmah.C")));

    // For Monte Carlo input
    if (isMC) {
        AliMCEventHandler *mcHandler = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mcHandler);
        mcHandler->SetPreReadMode(AliMCEventHandler::PreReadMode_t::kLmPreRead);
        // set below to kTRUE to read also track refs
        mcHandler->SetReadTR(kTRUE);
    }

    /// start analysis
    if (!mgr->InitAnalysis()) {
        return;
    }

    mgr->SetDebugLevel(3);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(kTRUE, 250);

    if (local) {

        /**********************/
        /*** LOCAL ANALYSIS ***/
        /**********************/

        myTask->EnableLocalMode();
        TChain *pChain = new TChain("esdTree");

        {
            if (sub == 702) {
                pChain->Add("/misc/alidata120/alice_u/schmah/TRD_self_tracking/Raw_data_pPb_2016/AliESDs.root");
            }
            if (sub == 800) {
                pChain->Add("/misc/alidata120/alice_u/schmah/TRD_self_tracking/AliESDs.root");
            }
            if (sub == 900) {
                pChain->Add("/misc/alidata131/alice_u/oschmidt/condor/sim_pp/chunk_982/AliESDs.root");
            }
            if (sub == 910) {
                pChain->Add("/misc/alidata100/alice_u/kroesen/tmp/AliESDs.root");
            }
            if (sub == 888) {
                pChain->Add("/misc/alidata141/alice_u/schmah/TRD/Calibration/Krypton/Data/raw/2018/AliESDs.root");
            }
            if (sub == 123) {
                pChain->Add("/misc/alidata141/alice_u/schmah/TRD/pp_sim/AliESDs.root");
            }
            if (sub == 555) {
                pChain->Add("AliESDs.root");
            }
        }

        mgr->StartAnalysis("local", pChain);
    } else {

        /**********************/
        /*** GRID ANALYSIS ***/
        /**********************/

        // set run list
        Int_t gridTest = 1;
        TString sGridMode(cGridMode);
        if (sGridMode == "full") {
            gridTest = 0;
        }
        if (sGridMode == "terminate") {
            gridTest = 2;
        }

        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();

        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

        // make sure your source files get copied to grid
        // alienHandler->SetAdditionalLibs(
        //     "Ali_AS_Event.h Ali_AS_EventLinkDef.h Ali_make_tracklets_from_digits.h Ali_make_tracklets_from_digits.cxx Ali_TRD_ST.h "
        //     "Ali_TRD_ST_LinkDef.h");
        alienHandler->SetAdditionalLibs(
            "Ali_AS_Event.h Ali_AS_Track.h Ali_Helix.h Ali_TRD_ST_LinkDef.h Ali_TRD_ST_Tracklets.h "
            "Ali_AS_EventLinkDef.h Ali_AS_Tracklet.h Ali_MC_particle.h Ali_TRD_ST_TOF_hit.h Ali_TRD_ST_V0.h "
            "Ali_AS_TRD_digit.h Ali_AS_offline_Tracklet.h Ali_TRD_ST_Event.h Ali_TRD_ST_TPC_Track.h "
            "Ali_make_tracklets_from_digits.h");
        alienHandler->SetAnalysisSource("Ali_make_tracklets_from_digits.cxx");

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)
        // alienHandler->SetAliPhysicsVersion("vAN-20200815_ROOT6-1");
        alienHandler->SetAliPhysicsVersion("vAN-20201103_ROOT6-1");
        // alienHandler->SetAliPhysicsVersion("vAN-20201203_ROOT6-1");
#else
        alienHandler->SetAliPhysicsVersion("vAN-20181112-1");
#endif

        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");

        // select the input data
        if (isMC && sub != 123) {
            alienHandler->SetGridDataDir("/alice/sim/2020/LHC20e3a/");
            alienHandler->SetDataPattern("*AOD/*AOD.root");
        } else {
            if (isMC && sub == 123) {
                alienHandler->SetGridDataDir("/alice/sim/2022/LHC22a3");
                alienHandler->SetGridWorkingDir(Form("%s/sub%d/", fname.Data(), sub));
                alienHandler->SetDataPattern("*/AliESDs.root");
                alienHandler->SetAnalysisMacro(Form("TaskTrackAna%d.C", sub));
                alienHandler->SetExecutable(Form("TaskTrackAna%d.sh", sub));
                alienHandler->SetJDLName(Form("TaskTrackAna%d.jdl", sub));
                alienHandler->SetRunPrefix("");

                Int_t runnumbers[] = {265309, 265332, 265334, 265336, 265338, 265339, 265342, 265343, 265344, 265377, 265378,
                                      265381, 265383, 265384, 265385, 265387, 265388, 265419, 265420, 265421, 265422, 265424,
                                      265425, 265426, 265427, 265435, 265499, 265500, 265501, 265521, 265525};  //  runs in total

                for (Int_t irun = 1; irun < 5; irun++) {
                    printf("%d %d", irun, runnumbers[irun]);
                    alienHandler->AddRunNumber(runnumbers[irun]);
                }
            }

            if (sub == 702) {
                alienHandler->SetGridDataDir("/alice/data/2016/LHC16q");
                alienHandler->SetGridWorkingDir(Form("%s/sub%d/", fname.Data(), sub));
                alienHandler->SetDataPattern("*pass1_CENT_wSDD/*/AliESDs.root");
                alienHandler->SetAnalysisMacro(Form("TaskTrackAna%d.C", sub));
                alienHandler->SetExecutable(Form("TaskTrackAna%d.sh", sub));
                alienHandler->SetJDLName(Form("TaskTrackAna%d.jdl", sub));
                alienHandler->SetRunPrefix("000");

                Int_t runnumbers[] = {265338, 265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424,
                                      265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378,
                                      265377, 265344, 265343, 265342, 265339, 265336, 265334, 265332, 265309};  // 31 runs in total

                for (Int_t irun = 0; irun < 1; irun++) {
                    printf("%d %d", irun, runnumbers[irun]);
                    alienHandler->AddRunNumber(runnumbers[irun]);
                }
            }
            if (sub == 7021) {
                // Full production for p-Pb, only one run: 265501
                alienHandler->SetGridDataDir("/alice/data/2016/LHC16q");
                alienHandler->SetGridWorkingDir(Form("%s/sub%d/", fname.Data(), sub));
                alienHandler->SetDataPattern("*pass2_TRD/*/AliESDs.root");
                alienHandler->SetAnalysisMacro(Form("TaskTrackAna%d.C", sub));
                alienHandler->SetExecutable(Form("TaskTrackAna%d.sh", sub));
                alienHandler->SetJDLName(Form("TaskTrackAna%d.jdl", sub));
                alienHandler->SetRunPrefix("000");

                Int_t runnumbers[] = {265501};  // 1 run in total

                for (Int_t irun = 0; irun < 1; irun++) {
                    printf("%d %d", irun, runnumbers[irun]);
                    alienHandler->AddRunNumber(runnumbers[irun]);
                }
            }
            if (sub == 800) {
                alienHandler->SetGridDataDir("/alice/data/2018/LHC18r");
                alienHandler->SetGridWorkingDir(Form("%s/sub%d/", fname.Data(), sub));
                alienHandler->SetDataPattern("*filter_trd_pass3/*/AliESDs.root");
                alienHandler->SetAnalysisMacro(Form("TaskTrackAna%d.C", sub));
                alienHandler->SetExecutable(Form("TaskTrackAna%d.sh", sub));
                alienHandler->SetJDLName(Form("TaskTrackAna%d.jdl", sub));
                alienHandler->SetRunPrefix("000");

                Int_t runnumbers[] = {296849, 296850, 297595};  // 3 runs in total

                for (Int_t irun = 0; irun < 1; irun++) {
                    printf("%d %d", irun, runnumbers[irun]);
                    alienHandler->AddRunNumber(runnumbers[irun]);
                }
            }
        }

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(50);  // 40
        // specify how many seconds your job may take
        alienHandler->SetTTL(86300);  // 39999

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);

        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(kFALSE)
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(useJDL);

        // define the output folders
        // alienHandler->SetGridOutputDir("myOutputDir_LHC18q_test1");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if (gridTest == 1) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        }
        if (gridTest == 0) {
            // else launch the full grid analysis
            alienHandler->SetRunMode("full");
            mgr->StartAnalysis("grid");
        }
        if (gridTest == 2) {
            // launch terminate to merge files
            alienHandler->SetRunMode("terminate");
            mgr->StartAnalysis("grid");
        }
    }
}
