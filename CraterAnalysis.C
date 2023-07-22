//
//      ********************
//      * CraterAnalysis.C *
//      ********************
//
//
//      The following program is developed for the evaluation of ASF-files from the SAMAICA-Software
//      to analyze crater distributions and their properties.
//
//      Programtypes:
//      -> CraterAnalysis: The absolute coordinates (centres) of craters are not controlled.
//      -> CraterAnalysisPosCheckEdge: The absolute coordinates of craters at the edge (overlap) are controlled.
//      -> CraterAnalysisPosCheckTotal: The absolute coordinates of craters over the whole area are controlled.
//
//      -------------------------------------------------------------------------------------------------------------
//      Input: Please note to use ASCII-Format (corresponds to data format 1 in Samaica).
//             Replace the commas as separators through blank spaces and delete the last unused line in the ASF-file.
//
//      Attention: If there are changes concerning the optimized values of the parameters lappx, lappy,
//                 xmin, xmax, ymin, ymax (acceptance window) please update these values in the CONSTRUCTOR!
//
//      Program execution: Open ROOT (5.28) and write the following commands:
//                         root [0] .L <path>\CraterAnalysis.C
//                         root [1] CraterAnalysis()
//      ____________________________________
//      Author: Johannes Knaute
//              April, 2011
//      ____________________________________
//


#include "TGButton.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TGTextEntry.h"
#include "TGDoubleSlider.h"
#include "TGLabel.h"
#include <TNtuple.h>

enum ETestCommandIdentifiers {
    HId1,
    HId2,
    HId3,
    HCId1,
    HCId2,
    HSId1
};

class TCraterAnalysis : public TGMainFrame {

private:

   //Elements For Defining Cuts:
    TGHorizontalFrame   *fFrame1_slider_x, *fFrame2_slider_x, *fFrame3_slider_x,//Slider1
                        *fFrame1_slider_y, *fFrame2_slider_y, *fFrame3_slider_y,//Slider2
                        *fFrame1_slider_b, *fFrame2_slider_b, *fFrame3_slider_b,//Slider3
                        *fFrame1_slider_e, *fFrame2_slider_e, *fFrame3_slider_e,//Slider4
                        *fFrame1_slider_area, *fFrame2_slider_area, *fFrame3_slider_area,//Slider5
                        *fFrame_button;
    TGLayoutHints       *fBuLy_slider_x, *fFrLy1_slider_x, *fFrLy2_slider_x, *fFrLy3_slider_x,//Slider1
                        *fBuLy_slider_y, *fFrLy1_slider_y, *fFrLy2_slider_y, *fFrLy3_slider_y,//Slider2
                        *fBuLy_slider_b, *fFrLy1_slider_b, *fFrLy2_slider_b, *fFrLy3_slider_b,//Slider3
                        *fBuLy_slider_e, *fFrLy1_slider_e, *fFrLy2_slider_e, *fFrLy3_slider_e,//Slider4
                        *fBuLy_slider_area, *fFrLy1_slider_area, *fFrLy2_slider_area, *fFrLy3_slider_area;//Slider5
    TGDoubleHSlider     *fslider_x,//Slider1
                        *fslider_y,//Slider2
                        *fSlider_b,//Slider3
                        *fSlider_e,//Slider4
                        *fSlider_area;//Slider5
    TGTextEntry         *fBoxEntry_xmin, *fBoxEntry_xmax,//Slider1
                        *fBoxEntry_ymin, *fBoxEntry_ymax,//Slider2
                        *fBoxEntry_bmin, *fBoxEntry_bmax,//Slider3
                        *fBoxEntry_emin, *fBoxEntry_emax,//Slider4
                        *fBoxEntry_areamin, *fBoxEntry_areamax,//Slider5
                        *fPath, *fFileName;//Edit-windows for analyzed Path and ASF-File
    TGTextBuffer        *fBoxBuffer_xmin, *fBoxBuffer_xmax,//Slider1
                        *fBoxBuffer_ymin, *fBoxBuffer_ymax,//Slider2
                        *fBoxBuffer_bmin, *fBoxBuffer_bmax,//Slider3
                        *fBoxBuffer_emin, *fBoxBuffer_emax,//Slider4
                        *fBoxBuffer_areamin, *fBoxBuffer_areamax;//Slider5
    TGLabel             *fLabel_x, *fLabel_y, *fLabel_b, *fLabel_e, *fLabel_area,//Labels for Explanations
                        *fLabel_path, *fLabel_filename;
    TGCheckButton       *fUnderground, *fPosCheck;//Check Button For Underground-Filtration
    TGVButtonGroup      *fButtonGroup;//Button Group For Position-Control-Method
    TGRadioButton       *fPosCheckMethod[2];//Radio Buttons For Choosing Position-Control-Method

   //Analysis - Definition Of Histograms:
    TH2F *hpos;
    TH2F *hpos_cut;
    TH1F *h_b_axis;
    TH1F *h_b_axis_cut;
    TH1F *hecc;
    TH1F *hecc_cut;
    TH1F *harea_ea;
    TH1F *harea_ca;
    TH1F *harea_ca_cut;
    TF1  *fgaus;

   //Cut-Variables:
    float fx_min;
    float fx_max;
    float fy_min;
    float fy_max;
    float fb_min;
    float fb_max;
    float fe_min;
    float fe_max;
    float farea_min;
    float farea_max;

   //Variable For Programtype:
    char program[200];

   //Variables For Path And Filename:
    char path[200];
    char filename[200];
    char datafile[200];

   //Number Of Retrievals:
    int counter;

   //Line-Number:
    int line;

   //Variables: Number Of Unfiltered And Filtered Craters:
    int Ncrater_unfiltered;
    int Ncrater_filtered;

   //Header 1 - Variables:
    int anzahl;
    float scalex,scaley;

   //Header 2 - Variables:
    int image_count_x,image_count_y;

   //Correlation Coefficients:
    float r_EA, r_EB, r_AB;

   //Scanarea:
    float A_scan;

   //Increments:
    float incx, incy;

   //Variables For The Overlap:
    float lappx, lappy;

   //Variables For The Acceptance Window:
    float xmin, xmax, ymin, ymax;

   //Definition Of A Tuple For Processing:
    TNtuple *tuple;

public:
    TCraterAnalysis();
    virtual ~TCraterAnalysis();
    void SetGroupEnabled(Bool_t);
    void CloseWindow();
    void DoText(const char *text);
    void DoSlider();
    void DoCanvas();
    void DoSave();

    ClassDef(TCraterAnalysis, 0)
};

//CONSTRUCTOR____________________________________________________________________________________________________________________________________________________________________
TCraterAnalysis::TCraterAnalysis() : TGMainFrame(gClient->GetRoot(), 100, 100)
{
    char buf[32];
    SetCleanup(kDeepCleanup);

   //Defining A Tuple For Processing The Data:
    tuple=new TNtuple("crater","crater","image_x:image_y:pos_u:pos_v:b:e:sphi:S2:ea:f_cb");

   //Analysis - Definition Of Histograms:
    hpos = new TH2F("hpos","Craterpositions",1400,-50,139950,1400,-50,139950);
    hpos_cut = new TH2F("hpos_cut","Craterpositions with Cuts",1400,-50,139950,1400,-50,139950);
    h_b_axis = new TH1F("h_b_axis","Semi Minor Axis",100,-0.02,3.98);
    h_b_axis_cut = new TH1F("h_b_axis_cut","Semi Minor Axis",100,-0.02,3.98);
    hecc = new TH1F("hecc","Eccentricity",220,-0.0025,1.0975);
    hecc_cut = new TH1F("hecc_cut","Eccentricity",220,-0.0025,1.0975);
    harea_ea = new TH1F("harea_ea","Enclosed Area",100,-0.2,39.8);
    harea_ca = new TH1F("harea_ca","Calculated Area",100,-0.2,39.8);
    harea_ca_cut = new TH1F("harea_ca_cut","Calculated Area",100,-0.2,39.8);
    fgaus=new TF1("fgaus","gaus",-0.2,40);

   //Defining Frames:
    fFrame1_slider_x = new TGHorizontalFrame(this, 0, 0, 0);//Slider1
    fFrame1_slider_x->Resize(200, 50);
    fFrame2_slider_x = new TGHorizontalFrame(this, 0, 0, 0);
    fFrame1_slider_y = new TGHorizontalFrame(this, 0, 0, 0);//Slider2
    fFrame1_slider_y->Resize(200, 50);
    fFrame2_slider_y = new TGHorizontalFrame(this, 0, 0, 0);
    fFrame1_slider_b = new TGHorizontalFrame(this, 0, 0, 0);//Slider3
    fFrame1_slider_b->Resize(200, 50);
    fFrame2_slider_b = new TGHorizontalFrame(this, 0, 0, 0);
    fFrame1_slider_e = new TGHorizontalFrame(this, 0, 0, 0);//Slider4
    fFrame1_slider_e->Resize(200, 50);
    fFrame2_slider_e = new TGHorizontalFrame(this, 0, 0, 0);
    fFrame1_slider_area = new TGHorizontalFrame(this, 0, 0, 0);//Slider5
    fFrame1_slider_area->Resize(200, 50);
    fFrame2_slider_area = new TGHorizontalFrame(this, 0, 0, 0);
    fFrame_button = new TGHorizontalFrame(this, 0, 0, 0);//Frame for button
    fFrame_button->Resize(110, 60);

   //Creating Sliders:
    fslider_x = new TGDoubleHSlider(fFrame2_slider_x, 410, kDoubleScaleBoth, HSId1,//Slider1
                                    kHorizontalFrame, GetDefaultFrameBackground(), kFALSE, kFALSE);
    fslider_y = new TGDoubleHSlider(fFrame2_slider_y, 410, kDoubleScaleBoth, HSId1,//Slider2
                                    kHorizontalFrame, GetDefaultFrameBackground(), kFALSE, kFALSE);
    fSlider_b = new TGDoubleHSlider(fFrame2_slider_b, 410, kDoubleScaleBoth, HSId1,//Slider3
                                    kHorizontalFrame, GetDefaultFrameBackground(), kFALSE, kFALSE);
    fSlider_e = new TGDoubleHSlider(fFrame2_slider_e, 410, kDoubleScaleBoth, HSId1,//Slider4
                                    kHorizontalFrame, GetDefaultFrameBackground(), kFALSE, kFALSE);
    fSlider_area = new TGDoubleHSlider(fFrame2_slider_area, 410, kDoubleScaleBoth, HSId1,//Slider5
                                       kHorizontalFrame, GetDefaultFrameBackground(), kFALSE, kFALSE);

   //Creating An OK-Button:
    TGTextButton *fOK_Button = new TGTextButton(fFrame_button,"&OK - Accept Cuts");
    fOK_Button->Connect("Clicked()","TCraterAnalysis",this,"DoCanvas()");
    fOK_Button->SetTextJustify(36);
    fOK_Button->SetMargins(0,0,0,0);
    fOK_Button->SetWrapLength(-1);
    fOK_Button->Resize(100,22);
    fOK_Button->MoveResize(5,5,100,22);
    ULong_t ucolor;
    gClient->GetColorByName("#cc0000",ucolor);
    fOK_Button->ChangeBackground(ucolor);
    fFrame_button->AddFrame(fOK_Button, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY));

   //Creating A Save-Button:
    TGTextButton *fSave_Button = new TGTextButton(fFrame_button,"&Save Results");
    fSave_Button->Connect("Clicked()","TCraterAnalysis",this,"DoSave()");
    fSave_Button->SetTextJustify(36);
    fSave_Button->SetMargins(0,0,0,0);
    fSave_Button->SetWrapLength(-1);
    fSave_Button->Resize(100,22);
    fSave_Button->MoveResize(5,30,100,22);
    ULong_t ucolor;
    gClient->GetColorByName("#6666cc",ucolor);
    fSave_Button->ChangeBackground(ucolor);
    fFrame_button->AddFrame(fSave_Button, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY));

   //Creating A Check-Button For Activation Of Underground Filtration:
    fUnderground = new TGCheckButton(this,"Apply Underground-Filtration");
    fUnderground->SetState(kButtonDown);
    AddFrame(fUnderground,  new TGLayoutHints(kLHintsLeft | kLHintsTop,150,0,10,0));

   //Creating A Check Button For Activation Of Position Control Via Button-Group:
    fPosCheck = new TGCheckButton(this,"Activate Position Check");
    AddFrame(fPosCheck,  new TGLayoutHints(kLHintsLeft | kLHintsTop,150,0,0,0));
    fButtonGroup = new TGVButtonGroup(this, "Position-Control-Method");
    fPosCheckMethod[0] = new TGRadioButton(fButtonGroup, new TGHotString("Edge-Objects"));
    fPosCheckMethod[1] = new TGRadioButton(fButtonGroup, new TGHotString("Total"));
    fButtonGroup->Show();
    AddFrame(fButtonGroup,  new TGLayoutHints(kLHintsLeft | kLHintsTop,159,0,0,0));
    fPosCheck->Connect("Toggled(Bool_t)", "TCraterAnalysis", this,"SetGroupEnabled(Bool_t)");
    fButtonGroup->SetRadioButtonExclusive(kTRUE);
    fPosCheckMethod[1]->SetOn();
    fButtonGroup->SetState(kFALSE);

   //Creating Labels:
    //Layout1:
    TGGC *fTextGC;
    const TGFont *font = gClient->GetFont("-*-times-bold-r-*-*-18-*-*-*-*-*-*-*");
    if (!font)
        font = gClient->GetResourcePool()->GetDefaultFont();
    FontStruct_t labelfont = font->GetFontStruct();
    GCValues_t   gval;
    gval.fMask = kGCBackground | kGCFont | kGCForeground;
    gval.fFont = font->GetFontHandle();
    gClient->GetColorByName("yellow", gval.fBackground);
    fTextGC = gClient->GetGC(&gval, kTRUE);
    ULong_t ycolor;
    gClient->GetColorByName("yellow", ycolor);

    //Layout2:
    TGGC *fTextGC2;
    const TGFont *font2 = gClient->GetFont("-*-times-bold-r-*-*-16-*-*-*-*-*-*-*");
    if (!font2)
        font2 = gClient->GetResourcePool()->GetDefaultFont();
    FontStruct_t labelfont2 = font2->GetFontStruct();
    GCValues_t   gval2;
    gval2.fMask = kGCBackground | kGCFont | kGCForeground;
    gval2.fFont = font2->GetFontHandle();
    gClient->GetColorByName("black", gval2.fBackground);
    fTextGC2 = gClient->GetGC(&gval2, kTRUE);
    ULong_t bcolor;
    gClient->GetColorByName("black", bcolor);

   //Characterization OF Sliders:
    //Connections:
    fslider_x->Connect("PositionChanged()", "TCraterAnalysis", this, "DoSlider()");//Slider1
    fslider_y->Connect("PositionChanged()", "TCraterAnalysis", this, "DoSlider()");//Slider2
    fSlider_b->Connect("PositionChanged()", "TCraterAnalysis", this, "DoSlider()");//Slider3
    fSlider_e->Connect("PositionChanged()", "TCraterAnalysis", this, "DoSlider()");//Slider4
    fSlider_area->Connect("PositionChanged()", "TCraterAnalysis", this, "DoSlider()");//Slider5

    //Range:
    fslider_x->SetRange(0,140000);//Slider1
    fslider_y->SetRange(0,140000);//Slider2
    fSlider_b->SetRange(0,4);//Slider3
    fSlider_e->SetRange(0,1);//Slider4
    fSlider_area->SetRange(0,40);//Slider5

   //Characterization OF Frames:
    fFrame2_slider_x->Resize(200, 25);//Slider1
    fFrame2_slider_y->Resize(200, 25);//Slider2
    fFrame2_slider_b->Resize(200, 25);//Slider3
    fFrame2_slider_e->Resize(200, 25);//Slider4
    fFrame2_slider_area->Resize(200, 25);//Slider5

    fFrame3_slider_x = new TGHorizontalFrame(this, 0, 0, 0);//Slider1
    fFrame3_slider_y = new TGHorizontalFrame(this, 0, 0, 0);//Slider2
    fFrame3_slider_b = new TGHorizontalFrame(this, 0, 0, 0);//Slider3
    fFrame3_slider_e = new TGHorizontalFrame(this, 0, 0, 0);//Slider4
    fFrame3_slider_area = new TGHorizontalFrame(this, 0, 0, 0);//Slider5

   //Characterization OF Boxes:
    //Initializations:
    fBoxEntry_xmin = new TGTextEntry(fFrame3_slider_x, fBoxBuffer_xmin = new TGTextBuffer(5), HId1);//Slider1
    fBoxEntry_xmax = new TGTextEntry(fFrame3_slider_x, fBoxBuffer_xmax = new TGTextBuffer(5), HId3);
    fBoxEntry_ymin = new TGTextEntry(fFrame3_slider_y, fBoxBuffer_ymin = new TGTextBuffer(5), HId1);//Slider2
    fBoxEntry_ymax = new TGTextEntry(fFrame3_slider_y, fBoxBuffer_ymax = new TGTextBuffer(5), HId3);
    fBoxEntry_bmin = new TGTextEntry(fFrame3_slider_b, fBoxBuffer_bmin = new TGTextBuffer(5), HId1);//Slider3
    fBoxEntry_bmax = new TGTextEntry(fFrame3_slider_b, fBoxBuffer_bmax = new TGTextBuffer(5), HId3);
    fBoxEntry_emin = new TGTextEntry(fFrame3_slider_e, fBoxBuffer_emin = new TGTextBuffer(5), HId1);//Slider4
    fBoxEntry_emax = new TGTextEntry(fFrame3_slider_e, fBoxBuffer_emax = new TGTextBuffer(5), HId3);
    fBoxEntry_areamin = new TGTextEntry(fFrame3_slider_area, fBoxBuffer_areamin = new TGTextBuffer(5), HId1);//Slider5
    fBoxEntry_areamax = new TGTextEntry(fFrame3_slider_area, fBoxBuffer_areamax = new TGTextBuffer(5), HId3);

    //Size:
    fBoxEntry_xmin->Resize(70,fBoxEntry_xmin->GetDefaultHeight());
    fBoxEntry_xmax->Resize(70,fBoxEntry_xmin->GetDefaultHeight());
    fBoxEntry_ymin->Resize(70,fBoxEntry_xmin->GetDefaultHeight());
    fBoxEntry_ymax->Resize(70,fBoxEntry_xmin->GetDefaultHeight());
    fBoxEntry_bmin->Resize(70,fBoxEntry_xmin->GetDefaultHeight());
    fBoxEntry_bmax->Resize(70,fBoxEntry_xmin->GetDefaultHeight());
    fBoxEntry_emin->Resize(70,fBoxEntry_xmin->GetDefaultHeight());
    fBoxEntry_emax->Resize(70,fBoxEntry_xmin->GetDefaultHeight());
    fBoxEntry_areamin->Resize(70,fBoxEntry_xmin->GetDefaultHeight());
    fBoxEntry_areamax->Resize(70,fBoxEntry_xmin->GetDefaultHeight());

    //Textinformation:
    fBoxEntry_xmin->SetToolTipText("x_min");//Slider1
    fBoxEntry_xmax->SetToolTipText("x_max");
    fBoxEntry_ymin->SetToolTipText("y_min");//Slider2
    fBoxEntry_ymax->SetToolTipText("y_max");
    fBoxEntry_bmin->SetToolTipText("b_min");//Slider3
    fBoxEntry_bmax->SetToolTipText("b_max");
    fBoxEntry_emin->SetToolTipText("e_min");//Slider4
    fBoxEntry_emax->SetToolTipText("e_max");
    fBoxEntry_areamin->SetToolTipText("area_min");//Slider5
    fBoxEntry_areamax->SetToolTipText("area_max");

    //Filling Buffers:
    fBoxBuffer_xmin->AddText(0, "0.0");//Slider1
    fBoxBuffer_xmax->AddText(0, "0.0");
    fBoxBuffer_ymin->AddText(0, "0.0");//Slider2
    fBoxBuffer_ymax->AddText(0, "0.0");
    fBoxBuffer_bmin->AddText(0, "0.0");//Slider3
    fBoxBuffer_bmax->AddText(0, "0.0");
    fBoxBuffer_emin->AddText(0, "0.0");//Slider4
    fBoxBuffer_emax->AddText(0, "0.0");
    fBoxBuffer_areamin->AddText(0, "0.0");//Slider5
    fBoxBuffer_areamax->AddText(0, "0.0");

    //Connections:
    fBoxEntry_xmin->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");//Slider1
    fBoxEntry_xmax->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");
    fBoxEntry_ymin->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");//Slider2
    fBoxEntry_ymax->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");
    fBoxEntry_bmin->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");//Slider3
    fBoxEntry_bmax->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");
    fBoxEntry_emin->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");//Slider4
    fBoxEntry_emax->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");
    fBoxEntry_areamin->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");//Slider5
    fBoxEntry_areamax->Connect("TextChanged(char*)", "TCraterAnalysis", this, "DoText(char*)");

   //Characterization OF Frames:
    fFrame3_slider_x->Resize(100, 25);//Slider1
    fFrame3_slider_y->Resize(100, 25);//Slider2
    fFrame3_slider_b->Resize(100, 25);//Slider3
    fFrame3_slider_e->Resize(100, 25);//Slider4
    fFrame3_slider_area->Resize(100, 25);//Slider5

   //Layout For Buttons: top align, equally expand horizontally
    fBuLy_slider_x = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 3, 3, 0, 0);//Slider1
    fBuLy_slider_y = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 3, 3, 0, 0);//Slider2
    fBuLy_slider_b = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 3, 3, 0, 0);//Slider3
    fBuLy_slider_e = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 3, 3, 0, 0);//Slider4
    fBuLy_slider_area = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 3, 3, 0, 0);//Slider5

   //Layout For The Frame: place at bottom, right aligned
    fFrLy1_slider_x = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 3, 3, 10, 0);//Slider1
    fFrLy2_slider_x = new TGLayoutHints(kLHintsTop | kLHintsLeft,    3, 3, 10, 0);
    fFrLy3_slider_x = new TGLayoutHints(kLHintsTop | kLHintsRight,   3, 3, 10, 0);
    fFrLy1_slider_y = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 3, 3, 10, 0);//Slider2
    fFrLy2_slider_y = new TGLayoutHints(kLHintsTop | kLHintsLeft,    3, 3, 10, 0);
    fFrLy3_slider_y = new TGLayoutHints(kLHintsTop | kLHintsRight,   3, 3, 10, 0);
    fFrLy1_slider_b = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 3, 3, 10, 0);//Slider3
    fFrLy2_slider_b = new TGLayoutHints(kLHintsTop | kLHintsLeft,    3, 3, 10, 0);
    fFrLy3_slider_b = new TGLayoutHints(kLHintsTop | kLHintsRight,   3, 3, 10, 0);
    fFrLy1_slider_e = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 3, 3, 10, 0);//Slider4
    fFrLy2_slider_e = new TGLayoutHints(kLHintsTop | kLHintsLeft,    3, 3, 10, 0);
    fFrLy3_slider_e = new TGLayoutHints(kLHintsTop | kLHintsRight,   3, 3, 10, 0);
    fFrLy1_slider_area = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 3, 3, 10, 5);//Slider5
    fFrLy2_slider_area = new TGLayoutHints(kLHintsTop | kLHintsLeft,    3, 3, 10, 5);
    fFrLy3_slider_area = new TGLayoutHints(kLHintsTop | kLHintsRight,   3, 3, 10, 5);

    fFrame2_slider_x->AddFrame(fslider_x, fBuLy_slider_x);//Slider1
    fFrame3_slider_x->AddFrame(fBoxEntry_xmin, fFrLy2_slider_x);
    fFrame3_slider_x->AddFrame(fBoxEntry_xmax, fFrLy3_slider_x);
    fFrame2_slider_y->AddFrame(fslider_y, fBuLy_slider_y);//Slider2
    fFrame3_slider_y->AddFrame(fBoxEntry_ymin, fFrLy2_slider_y);
    fFrame3_slider_y->AddFrame(fBoxEntry_ymax, fFrLy3_slider_y);
    fFrame2_slider_b->AddFrame(fSlider_b, fBuLy_slider_b);//Slider3
    fFrame3_slider_b->AddFrame(fBoxEntry_bmin, fFrLy2_slider_b);
    fFrame3_slider_b->AddFrame(fBoxEntry_bmax, fFrLy3_slider_b);
    fFrame2_slider_e->AddFrame(fSlider_e, fBuLy_slider_e);//Slider4
    fFrame3_slider_e->AddFrame(fBoxEntry_emin, fFrLy2_slider_e);
    fFrame3_slider_e->AddFrame(fBoxEntry_emax, fFrLy3_slider_e);
    fFrame2_slider_area->AddFrame(fSlider_area, fBuLy_slider_area);//Slider5
    fFrame3_slider_area->AddFrame(fBoxEntry_areamin, fFrLy2_slider_area);
    fFrame3_slider_area->AddFrame(fBoxEntry_areamax, fFrLy3_slider_area);

   //Drawing Elements: Labels + Sliders:
    fLabel_x = new TGLabel(this, "Position x [µm]", fTextGC->GetGC(),labelfont);//Label1
    AddFrame(fLabel_x,  new TGLayoutHints(kLHintsCenterX, 0, 0, 0, 0));
    fLabel_x->SetTextColor(ycolor);
    AddFrame(fFrame1_slider_x, fBuLy_slider_x);//Slider1
    AddFrame(fFrame2_slider_x, fBuLy_slider_x);
    AddFrame(fFrame3_slider_x, fBuLy_slider_x);

    fLabel_y = new TGLabel(this, "Position y [µm]", fTextGC->GetGC(),labelfont);//Label2
    AddFrame(fLabel_y,  new TGLayoutHints(kLHintsCenterX, 0, 0, 15, 0));
    fLabel_y->SetTextColor(ycolor);
    AddFrame(fFrame1_slider_y, fBuLy_slider_y);//Slider2
    AddFrame(fFrame2_slider_y, fBuLy_slider_y);
    AddFrame(fFrame3_slider_y, fBuLy_slider_y);

    fLabel_b = new TGLabel(this, "Semi Minor Axis [µm]", fTextGC->GetGC(),labelfont);//Label3
    AddFrame(fLabel_b,  new TGLayoutHints(kLHintsCenterX, 0, 0, 15, 0));
    fLabel_b->SetTextColor(ycolor);
    AddFrame(fFrame1_slider_b, fBuLy_slider_b);//Slider3
    AddFrame(fFrame2_slider_b, fBuLy_slider_b);
    AddFrame(fFrame3_slider_b, fBuLy_slider_b);

    fLabel_e = new TGLabel(this, "Eccentricity [1]", fTextGC->GetGC(),labelfont);//Label4
    AddFrame(fLabel_e,  new TGLayoutHints(kLHintsCenterX, 0, 0, 15, 0));
    fLabel_e->SetTextColor(ycolor);
    AddFrame(fFrame1_slider_e, fBuLy_slider_e);//Slider4
    AddFrame(fFrame2_slider_e, fBuLy_slider_e);
    AddFrame(fFrame3_slider_e, fBuLy_slider_e);

    fLabel_area = new TGLabel(this, "Calculated Area [µm²]", fTextGC->GetGC(),labelfont);//Label5
    AddFrame(fLabel_area,  new TGLayoutHints(kLHintsCenterX, 0, 0, 15, 0));
    fLabel_area->SetTextColor(ycolor);
    AddFrame(fFrame1_slider_area, fBuLy_slider_area);//Slider5
    AddFrame(fFrame2_slider_area, fBuLy_slider_area);
    AddFrame(fFrame3_slider_area, fBuLy_slider_area);

   //Creating Edit-Windows For Path And File-Name With Accompanying Labels:
    //Path-Entry:
    fLabel_path = new TGLabel(this, " Path of File for Analysis: [Form:   C:/Users/.../.../]", fTextGC2->GetGC(),labelfont2);
    AddFrame(fLabel_path,  new TGLayoutHints(kLHintsNormal, 0, 0, 15, 0));
    fLabel_path->SetTextColor(bcolor);
    fPath = new TGTextEntry(this);
    AddFrame(fPath, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));
    fPath->SetText("path");

    //Filename-Entry:
    fLabel_filename = new TGLabel(this, " Name of ASF-File: [Form:   abc]", fTextGC2->GetGC(),labelfont2);
    AddFrame(fLabel_filename,  new TGLayoutHints(kLHintsNormal, 0, 0, 15, 0));
    fLabel_filename->SetTextColor(bcolor);
    fFileName = new TGTextEntry(this);
    AddFrame(fFileName, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));
    fFileName->SetText("filename");

   //Set Main Frame Name, Map Sub Windows (Buttons), Initialize Layout
    //Algorithm Via Resize() And Map Main Frame
    SetWindowName("DEFINE CUTS");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

   //Defining Start-Values OF Sliders:
    fslider_x->SetPosition(0,140000);//Slider1
    fslider_y->SetPosition(0,140000);//Slider2
    fSlider_b->SetPosition(0,4);//Slider3
    fSlider_e->SetPosition(0,1);//Slider4
    fSlider_area->SetPosition(0,40);//Slider5

   //Writing Values into Boxes:
    sprintf(buf, "%.0f", fslider_x->GetMinPosition());//Slider1
    fBoxBuffer_xmin->Clear();
    fBoxBuffer_xmin->AddText(0, buf);
    sprintf(buf, "%.0f", fslider_x->GetMaxPosition());
    fBoxBuffer_xmax->Clear();
    fBoxBuffer_xmax->AddText(0, buf);
    sprintf(buf, "%.0f", fslider_y->GetMinPosition());//Slider2
    fBoxBuffer_ymin->Clear();
    fBoxBuffer_ymin->AddText(0, buf);
    sprintf(buf, "%.0f", fslider_y->GetMaxPosition());
    fBoxBuffer_ymax->Clear();
    fBoxBuffer_ymax->AddText(0, buf);
    sprintf(buf, "%.3f", fSlider_b->GetMinPosition());//Slider3
    fBoxBuffer_bmin->Clear();
    fBoxBuffer_bmin->AddText(0, buf);
    sprintf(buf, "%.3f", fSlider_b->GetMaxPosition());
    fBoxBuffer_bmax->Clear();
    fBoxBuffer_bmax->AddText(0, buf);
    sprintf(buf, "%.3f", fSlider_e->GetMinPosition());//Slider4
    fBoxBuffer_emin->Clear();
    fBoxBuffer_emin->AddText(0, buf);
    sprintf(buf, "%.3f", fSlider_e->GetMaxPosition());
    fBoxBuffer_emax->Clear();
    fBoxBuffer_emax->AddText(0, buf);
    sprintf(buf, "%.3f", fSlider_area->GetMinPosition());//Slider5
    fBoxBuffer_areamin->Clear();
    fBoxBuffer_areamin->AddText(0, buf);
    sprintf(buf, "%.3f", fSlider_area->GetMaxPosition());
    fBoxBuffer_areamax->Clear();
    fBoxBuffer_areamax->AddText(0, buf);

   //Defining Cut-Values:
    fx_min=fslider_x->GetMinPosition();
    fx_max=fslider_x->GetMaxPosition();
    fy_min=fslider_y->GetMinPosition();
    fy_max=fslider_y->GetMaxPosition();
    fb_min=fSlider_b->GetMinPosition();
    fb_max=fSlider_b->GetMaxPosition();
    fe_min=fSlider_e->GetMinPosition();
    fe_max=fSlider_e->GetMaxPosition();
    farea_min=fSlider_area->GetMinPosition();
    farea_max=fSlider_area->GetMaxPosition();

   //Number of retrievals:
    counter=0;

   //Values Of The Overlap:
    lappx=60;//written down in the parameter file!
    lappy=60;

   //Values Of The Acceptance Window:
    xmin=-482;//written down in the parameter file!
    xmax=482;
    ymin=-482;
    ymax=482;
}

//DESTRUCTOR_____________________________________________________________________________________________________________________________________________________________________
TCraterAnalysis::~TCraterAnalysis()
{
    Cleanup();
}

//Enable The Button Group__FUNCTION_SET_GROUP_ENABLED____________________________________________________________________________________________________________________________
void TCraterAnalysis::SetGroupEnabled(Bool_t on)
{
    fButtonGroup->SetState(on);
}

//Closing The Window__FUNCTION_CLOSE_WINDOW______________________________________________________________________________________________________________________________________
void TCraterAnalysis::CloseWindow()
{
    delete this;//Called when window is closed via the window manager.
}

//Realization Of The Text-Entries__FUNCTION_DO_TEXT______________________________________________________________________________________________________________________________
void TCraterAnalysis::DoText(const char * /*text*/)
{
   //Handle Text Entry Widgets:
    //Minimum-Value:
    fslider_x->SetPosition(atof(fBoxBuffer_xmin->GetString()), fslider_x->GetMaxPosition());//Slider1
    fslider_y->SetPosition(atof(fBoxBuffer_ymin->GetString()), fslider_y->GetMaxPosition());//Slider2
    fSlider_b->SetPosition(atof(fBoxBuffer_bmin->GetString()), fSlider_b->GetMaxPosition());//Slider3
    fSlider_e->SetPosition(atof(fBoxBuffer_emin->GetString()), fSlider_e->GetMaxPosition());//Slider4
    fSlider_area->SetPosition(atof(fBoxBuffer_areamin->GetString()), fSlider_area->GetMaxPosition());//Slider5
    fx_max=atof(fBoxBuffer_xmax->GetString());//Slider1
    fx_min=atof(fBoxBuffer_xmin->GetString());
    fy_max=atof(fBoxBuffer_ymax->GetString());//Slider2
    fy_min=atof(fBoxBuffer_ymin->GetString());
    fb_max=atof(fBoxBuffer_bmax->GetString());//Slider3
    fb_min=atof(fBoxBuffer_bmin->GetString());
    fe_max=atof(fBoxBuffer_emax->GetString());//Slider4
    fe_min=atof(fBoxBuffer_emin->GetString());
    farea_max=atof(fBoxBuffer_areamax->GetString());//Slider5
    farea_min=atof(fBoxBuffer_areamin->GetString());

    //Maximum-Value:
    fslider_x->SetPosition(fslider_x->GetMinPosition(), atof(fBoxBuffer_xmax->GetString()));//Slider1
    fslider_y->SetPosition(fslider_y->GetMinPosition(), atof(fBoxBuffer_ymax->GetString()));//Slider2
    fSlider_b->SetPosition(fSlider_b->GetMinPosition(), atof(fBoxBuffer_bmax->GetString()));//Slider3
    fSlider_e->SetPosition(fSlider_e->GetMinPosition(), atof(fBoxBuffer_emax->GetString()));//Slider4
    fSlider_area->SetPosition(fSlider_area->GetMinPosition(), atof(fBoxBuffer_areamax->GetString()));//Slider5
    fx_max=atof(fBoxBuffer_xmax->GetString());//Slider1
    fx_min=atof(fBoxBuffer_xmin->GetString());
    fy_max=atof(fBoxBuffer_ymax->GetString());//Slider2
    fy_min=atof(fBoxBuffer_ymin->GetString());
    fb_max=atof(fBoxBuffer_bmax->GetString());//Slider3
    fb_min=atof(fBoxBuffer_bmin->GetString());
    fe_max=atof(fBoxBuffer_emax->GetString());//Slider4
    fe_min=atof(fBoxBuffer_emin->GetString());
    farea_max=atof(fBoxBuffer_areamax->GetString());//Slider5
    farea_min=atof(fBoxBuffer_areamin->GetString());

    cout<<"DoText"<<endl;
}

//Realization Of The Slider Movements__FUNCTION_DO_SLIDER________________________________________________________________________________________________________________________
void TCraterAnalysis::DoSlider()
{
    char buf[32];

   //Reading Values:
    fx_min=fslider_x->GetMinPosition();//Slider1
    fx_max=fslider_x->GetMaxPosition();
    fy_min=fslider_y->GetMinPosition();//Slider2
    fy_max=fslider_y->GetMaxPosition();
    fb_min=fSlider_b->GetMinPosition();//Slider3
    fb_max=fSlider_b->GetMaxPosition();
    fe_min=fSlider_e->GetMinPosition();//Slider4
    fe_max=fSlider_e->GetMaxPosition();
    farea_min=fSlider_area->GetMinPosition();//Slider5
    farea_max=fSlider_area->GetMaxPosition();

   //Realizations Of The Movements:
    sprintf(buf, "%.0f", fslider_x->GetMinPosition());//Slider1
    fBoxBuffer_xmin->Clear();
    fBoxBuffer_xmin->AddText(0, buf);
    fBoxEntry_xmin->SetCursorPosition(fBoxEntry_xmin->GetCursorPosition());
    fBoxEntry_xmin->Deselect();
    gClient->NeedRedraw(fBoxEntry_xmin);
    sprintf(buf, "%.0f", fslider_x->GetMaxPosition());
    fBoxBuffer_xmax->Clear();
    fBoxBuffer_xmax->AddText(0, buf);
    fBoxEntry_xmax->SetCursorPosition(fBoxEntry_xmax->GetCursorPosition());
    fBoxEntry_xmax->Deselect();
    gClient->NeedRedraw(fBoxEntry_xmax);

    sprintf(buf, "%.0f", fslider_y->GetMinPosition());//Slider2
    fBoxBuffer_ymin->Clear();
    fBoxBuffer_ymin->AddText(0, buf);
    fBoxEntry_ymin->SetCursorPosition(fBoxEntry_ymin->GetCursorPosition());
    fBoxEntry_ymin->Deselect();
    gClient->NeedRedraw(fBoxEntry_ymin);
    sprintf(buf, "%.0f", fslider_y->GetMaxPosition());
    fBoxBuffer_ymax->Clear();
    fBoxBuffer_ymax->AddText(0, buf);
    fBoxEntry_ymax->SetCursorPosition(fBoxEntry_ymax->GetCursorPosition());
    fBoxEntry_ymax->Deselect();
    gClient->NeedRedraw(fBoxEntry_ymax);

    sprintf(buf, "%.3f", fSlider_b->GetMinPosition());//Slider3
    fBoxBuffer_bmin->Clear();
    fBoxBuffer_bmin->AddText(0, buf);
    fBoxEntry_bmin->SetCursorPosition(fBoxEntry_bmin->GetCursorPosition());
    fBoxEntry_bmin->Deselect();
    gClient->NeedRedraw(fBoxEntry_bmin);
    sprintf(buf, "%.3f", fSlider_b->GetMaxPosition());
    fBoxBuffer_bmax->Clear();
    fBoxBuffer_bmax->AddText(0, buf);
    fBoxEntry_bmax->SetCursorPosition(fBoxEntry_bmax->GetCursorPosition());
    fBoxEntry_bmax->Deselect();
    gClient->NeedRedraw(fBoxEntry_bmax);

    sprintf(buf, "%.3f", fSlider_e->GetMinPosition());//Slider4
    fBoxBuffer_emin->Clear();
    fBoxBuffer_emin->AddText(0, buf);
    fBoxEntry_emin->SetCursorPosition(fBoxEntry_emin->GetCursorPosition());
    fBoxEntry_emin->Deselect();
    gClient->NeedRedraw(fBoxEntry_emin);
    sprintf(buf, "%.3f", fSlider_e->GetMaxPosition());
    fBoxBuffer_emax->Clear();
    fBoxBuffer_emax->AddText(0, buf);
    fBoxEntry_emax->SetCursorPosition(fBoxEntry_emax->GetCursorPosition());
    fBoxEntry_emax->Deselect();
    gClient->NeedRedraw(fBoxEntry_emax);

    sprintf(buf, "%.3f", fSlider_area->GetMinPosition());//Slider5
    fBoxBuffer_areamin->Clear();
    fBoxBuffer_areamin->AddText(0, buf);
    fBoxEntry_areamin->SetCursorPosition(fBoxEntry_areamin->GetCursorPosition());
    fBoxEntry_areamin->Deselect();
    gClient->NeedRedraw(fBoxEntry_areamin);
    sprintf(buf, "%.3f", fSlider_area->GetMaxPosition());
    fBoxBuffer_areamax->Clear();
    fBoxBuffer_areamax->AddText(0, buf);
    fBoxEntry_areamax->SetCursorPosition(fBoxEntry_areamax->GetCursorPosition());
    fBoxEntry_areamax->Deselect();
    gClient->NeedRedraw(fBoxEntry_areamax);

    cout<<"DoSlider"<<endl;
}

//Creating Canvas And Processing Calculations__FUNCTION_DO_CANVAS________________________________________________________________________________________________________________
void TCraterAnalysis::DoCanvas()
{

//GLOBAL AREA////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    counter=counter+1;//Number of retrievals
    cout<< counter <<endl;
    line=0;

    ifstream asffile;

   //Chart Options:
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(1);
    gStyle->SetPalette(8);

   //Variables: Number Of Unfiltered And Filtered Craters:
    Ncrater_unfiltered=0;
    Ncrater_filtered=0;

   //Variables For Underground-Calculation:
    int K_v=0;
    float D_kv=0;
    float A_v=0;
    float A_a=0;
    int K_ua=0;

   //Variables For Filterlevels:
    int PositionFilter=0;
    int AxisFilter=0;
    int EccentricityFilter=0;
    int AreaFilter=0;

   //Correlation Coefficients:
    r_EA=0;//Eccentricity--Area
    r_EB=0;//Eccentricity--Semi Minor Axis
    r_AB=0;//Area--Semi Minor Axis
    float x_E=0;
    float x_A=0;
    float y_A=0;
    float y_B=0;
    float xy_EA=0;
    float xy_EB=0;
    float xy_AB=0;
    float x2_E=0;
    float x2_A=0;
    float y2_A=0;
    float y2_B=0;

   //Header 1 - Variables:
    anzahl=0;
    scalex=0;
    scaley=0;
    int ibright, icontrast;
    char dummy[200];

   //Header 2 - Variables:
    image_count_x=0;
    image_count_y=0;
    incx=0;
    incy=0;
    int mark1x,mark1y,mark2x,mark2y,startx,starty;

   //Array For Absolute Coordinates:
    float positionx[20000000];
    float positiony[20000000];

//BODY///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //Delete Old Histograms:
    if(counter>1){
       hpos->Reset();
       hpos_cut->Reset();
       h_b_axis->Reset();
       h_b_axis_cut->Reset();
       hecc->Reset();
       hecc_cut->Reset();
       harea_ea->Reset();
       harea_ca->Reset();
       harea_ca_cut->Reset();
    }

   //Output-Variables:
    int image_x,image_y;
    float pos_u,pos_v,b,e,sphi,S2,ea,f_cb;

   //File For Analysis:
    sprintf(path,fPath->GetText());
    sprintf(filename,fFileName->GetText());
    sprintf (datafile, "%s%s%s", path, filename, ".ASF");
    asffile.open(datafile);
    if (asffile.good()){

        //Header 1 - Reading Variables:
        asffile>>anzahl;
        asffile>>scalex>>scaley;
        asffile>>ibright>>icontrast;
        asffile>>dummy>>dummy>>dummy>>dummy;

        //Header 2 - Reading Variables:
        asffile>>image_count_x>>image_count_y>>mark1x>>mark1y>>mark2x>>mark2y>>incx>>incy>>startx>>starty;
    }
    else {cout<<"File does not exist"<<endl;}

   //Redefining Increments:
    incx=(1024-lappx)/fabs(scalex);
    incy=(1024-lappy)/scaley;

//CALCULATIONS///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    while(asffile.good())
    {//CALCULATION_LOOP//START///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //Line-Number:
        line=line+1;
        if ((fPosCheck->IsOn() && fPosCheckMethod[0]->IsOn()) || (fPosCheck->IsOn() && fPosCheckMethod[1]->IsOn()))
            cout<<line<<endl;

        //Reading Values:
        asffile>>image_x>>image_y>>pos_u>>pos_v>>b>>e>>sphi>>S2>>ea>>f_cb;
        tuple->Fill(image_x,image_y,pos_u,pos_v,b,e,sphi,S2,ea,f_cb);

        //Calculation Of Absolute Coordinates X,Y [µm]:
        float X,Y;
        X=(startx + image_x*incx) + (pos_u/scalex);
        Y=(starty + image_y*incy) + (pos_v/scaley);

        //Exclusion Of Same Objects:
        //CraterAnalysis:
        int i;
        float valuex, valuey;
        bool posx=false;
        bool posy=false;

        //CraterAnalysisPosCheckEdge:
        if (fPosCheck->IsOn() && fPosCheckMethod[0]->IsOn())
        {
            positionx[line-1]=X;
            positiony[line-1]=Y;
            if (pos_u<=-(512-lappx) || pos_u>=(512-lappx) || pos_v<=-(512-lappy) || pos_v>=(512-lappy))//Control only if it is an object at the edge of the video frame (time control)
            {
                for (i=0;i<(line-1);i++){
                    valuex=positionx[i];
                    if (X>=(valuex-0.5) && X<=(valuex+0.5)){
                        posx=true;
                        valuey=positiony[i];
                        if (Y>=(valuey-0.5) && Y<=(valuey+0.5)){
                            posy=true;
                            break;
                        }
                    }
                }
            }
        }

        //CraterAnalysisPosCheckTotal:
        if (fPosCheck->IsOn() && fPosCheckMethod[1]->IsOn())
        {
            positionx[line-1]=X;
            positiony[line-1]=Y;
            for (i=0;i<(line-1);i++){
                valuex=positionx[i];
                if (X>=(valuex-0.5) && X<=(valuex+0.5)){
                    posx=true;
                    valuey=positiony[i];
                    if (Y>=(valuey-0.5) && Y<=(valuey+0.5)){
                        posy=true;
                        break;
                    }
                }
            }
        }

        if (!(posx==true && posy==true))
        {
            Ncrater_unfiltered=Ncrater_unfiltered+1;

            //Calculation Of Major Axis And Converting Axis Into µm-Length:
            const double Pi = 3.1415926535897932384626433832795;
            float phi, q, r, s, t, a, a_new, b_new;
            a=b/e;
            phi=asin(sphi)*180.0/Pi;
            q=fabs(scaley/scalex);
            r=q*( (((cos(phi*Pi/180.0))*(cos(phi*Pi/180.0)))/pow(a,2)) + (((sin(phi*Pi/180.0))*(sin(phi*Pi/180.0)))/pow(b,2)) );
            s=(1/q) * ( (((sin(phi*Pi/180.0))*(sin(phi*Pi/180.0)))/pow(a,2)) + (((cos(phi*Pi/180.0))*(cos(phi*Pi/180.0)))/pow(b,2)) );
            t=( (1/pow(a,2)) - (1/pow(b,2)) ) * sin(phi*Pi/180.0) * cos(phi*Pi/180.0);
            a_new=sqrt(2/(r + s - sqrt( ((r-s)*(r-s))+4*t*t ) ) )   *sqrt(q)/fabs(scalex);
            b_new=sqrt(2/(r + s + sqrt( ((r-s)*(r-s))+4*t*t ) ) )   *sqrt(q)/fabs(scalex);

            //Conversion Factor For µm^2:
            float cfm;
            cfm=1/(fabs(scalex)*scaley);

            //Filling Standard-Histograms:
            float ea_new;
            hpos->Fill(X,Y);
            h_b_axis->Fill(b_new);
            hecc->Fill(e);
            ea_new=ea*cfm;
            harea_ea->Fill(ea_new);

            //Calculation Of Ellipse-Area And Histogram Of Calculated Area:
            float ca,ca_cut;
            ca=Pi*a_new*b_new;
            harea_ca->Fill(ca);
            ca_cut=ca;

            //Count For Filterlevels:
            if (X>=fx_min && X<=fx_max && Y>=fy_min && Y<=fy_max)
                PositionFilter=PositionFilter+1;
            if (b_new>=fb_min && b_new<=fb_max)
                AxisFilter=AxisFilter+1;
            if (e>=fe_min && e<=fe_max)
                EccentricityFilter=EccentricityFilter+1;
            if (ca>=farea_min && ca<=farea_max)
                AreaFilter=AreaFilter+1;

            //Count For Underground:
            if (X<fx_min || X>fx_max || Y<fy_min || Y>fy_max){
                if (b_new>=fb_min && b_new<=fb_max && e>=fe_min && e<=fe_max && ca>=farea_min && ca<=farea_max){
                    K_v=K_v+1;
                }
            }

            //Calculations For Correlation Coefficients:
            //Calculation Of r_EA:
            if (e<fe_min || e>fe_max){
                x_E=x_E+1;
                x2_E=x2_E+1;
            }
            if (ca<farea_min || ca>farea_max){
                y_A=y_A+1;
                y2_A=y2_A+1;
            }
            if ((e<fe_min || e>fe_max) && (ca<farea_min || ca>farea_max)){
                xy_EA=xy_EA+1;
            }

            //Additional Calculations For r_EB:
            if (b_new<fb_min || b_new>fb_max){
                y_B=y_B+1;
                y2_B=y2_B+1;
            }
            if ((e<fe_min || e>fe_max) && (b_new<fb_min || b_new>fb_max)){
                xy_EB=xy_EB+1;
            }

            //Additional Calculations For r_AB:
            if (ca<farea_min || ca>farea_max){
                x_A=x_A+1;
                x2_A=x2_A+1;
            }
            if ((ca<farea_min || ca>farea_max) && (b_new<fb_min || b_new>fb_max)){
                xy_AB=xy_AB+1;
            }

            //Creating Cut-Histograms:
            if (ca>=farea_min && ca<=farea_max && e>=fe_min && e<=fe_max && X>=fx_min && X<=fx_max && Y>=fy_min && Y<=fy_max && b_new>=fb_min && b_new<=fb_max){
                Ncrater_filtered=Ncrater_filtered+1;
                hecc_cut->Fill(e);
                hpos_cut->Fill(X,Y);
                h_b_axis_cut->Fill(b_new);
                harea_ca_cut->Fill(ca_cut);
            }
        }
    }//CALCULATION_LOOP//END/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    asffile.close();
//OTHER_CALCULATIONS/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //Calculation Of Correlation Coefficients:
    float n;
    n=Ncrater_unfiltered;
    r_EA=((xy_EA/n)-((x_E/n)*(y_A/n))) / ( (sqrt((x2_E/n)-pow((x_E/n),2)))  *  (sqrt((y2_A/n)-pow((y_A/n),2))) );
    r_EB=((xy_EB/n)-((x_E/n)*(y_B/n))) / ( (sqrt((x2_E/n)-pow((x_E/n),2)))  *  (sqrt((y2_B/n)-pow((y_B/n),2))) );
    r_AB=((xy_AB/n)-((x_A/n)*(y_B/n))) / ( (sqrt((x2_A/n)-pow((x_A/n),2)))  *  (sqrt((y2_B/n)-pow((y_B/n),2))) );

   //Calculation Of Scanned Area:
    float L1,L2;
    L1=image_count_x*incx + (xmax-xmin)/fabs(scalex);
    L2=image_count_y*incy + (ymax-ymin)/scaley;
    A_scan=L1*L2;//[µm²]

   //Calculation OF Filtered Craternumber Without Underground In Case Of Activation:
    if (fUnderground->IsOn())
    {
        A_a=(fx_max-fx_min)*(fy_max-fy_min);
        A_v=A_scan-A_a;
        D_kv=K_v/A_v;
        K_ua=D_kv*A_a;
        Ncrater_filtered=Ncrater_filtered-fabs(K_ua);
    }

//OUTPUT//START//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //Creating A Canvas:
    TCanvas *coverview = new TCanvas("coverview","Overview: CraterAnalysis",1300,700);
    coverview->Divide(4,2);

   //Filling Canvas:
    coverview->cd(1);//first window in canvas
        hpos->GetXaxis()->SetTitle("Position x [µm]");
        hpos->GetYaxis()->SetTitle("Position y [µm]");
        gStyle->SetPalette(1);
        hpos->Draw("colz");
    coverview->cd(2);//second window in canvas
        h_b_axis_cut->GetXaxis()->SetTitle("b [µm]");
        h_b_axis_cut->SetLineColor(2);
        h_b_axis_cut->SetLineWidth(2);
        h_b_axis_cut->Draw();
        h_b_axis->SetLineColor(1);
        h_b_axis->SetLineStyle(3);
        h_b_axis->Draw("same");
    coverview->cd(3);//third window in canvas
        harea_ea->GetXaxis()->SetTitle("enclosed area [µm^{2}]");
        harea_ea->Draw();
    coverview->cd(4);//fourth window in canvas
        TLatex Tl;
        Tl.SetTextAlign(12);
        Tl.SetTextSize(0.04);
        if (!(fPosCheck->IsOn()))
            Tl.DrawLatex(.05,.95,"Parameters (CraterAnalysis):");
        if (fPosCheck->IsOn() && fPosCheckMethod[0]->IsOn())
            Tl.DrawLatex(.05,.95,"Parameters (CraterAnalysisPosCheckEdge):");
        if (fPosCheck->IsOn() && fPosCheckMethod[1]->IsOn())
            Tl.DrawLatex(.05,.95,"Parameters (CraterAnalysisPosCheckTotal):");
        Tl.DrawLatex(0.05,0.857,"1)   Number Of Unfiltered Craters = ");
            sprintf(dummy,"%4.0f",Ncrater_unfiltered);
            TText *tparameter1=new TText(.7,.85,dummy);
            tparameter1->Draw();
        Tl.DrawLatex(0.05,0.757,"2)   Number Of Filtered Craters = ");
            sprintf(dummy,"%4.0f",Ncrater_filtered);
            TText *tparameter2=new TText(.7,.75,dummy);
            tparameter2->Draw();
        Tl.DrawLatex(0.05,0.65,"Cuts:");
        Tl.DrawLatex(0.05,0.557,"1)   x_min, x_max: ");
            sprintf(dummy,"%4.0f",fx_min);
            TText *tparameter3=new TText(.5,.55,dummy);
            tparameter3->Draw();
            sprintf(dummy,"%4.0f",fx_max);
            TText *tparameter4=new TText(.7,.55,dummy);
            tparameter4->Draw();
        Tl.DrawLatex(0.05,0.457,"2)   y_min, y_max: ");
            sprintf(dummy,"%4.0f",fy_min);
            TText *tparameter5=new TText(.5,.45,dummy);
            tparameter5->Draw();
            sprintf(dummy,"%4.0f",fy_max);
            TText *tparameter6=new TText(.7,.45,dummy);
            tparameter6->Draw();
        Tl.DrawLatex(0.05,0.357,"3)   b_min, b_max: ");
            sprintf(dummy,"%4.2f",fb_min);
            TText *tparameter7=new TText(.5,.35,dummy);
            tparameter7->Draw();
            sprintf(dummy,"%4.2f",fb_max);
            TText *tparameter8=new TText(.7,.35,dummy);
            tparameter8->Draw();
        Tl.DrawLatex(0.05,0.257,"4)   e_min, e_max: ");
            sprintf(dummy,"%4.2f",fe_min);
            TText *tparameter9=new TText(.5,.25,dummy);
            tparameter9->Draw();
            sprintf(dummy,"%4.2f",fe_max);
            TText *tparameter10=new TText(.7,.25,dummy);
            tparameter10->Draw();
        Tl.DrawLatex(0.05,0.157,"5)   area_min, area_max: ");
            sprintf(dummy,"%4.1f",farea_min);
            TText *tparameter11=new TText(.5,.15,dummy);
            tparameter11->Draw();
            sprintf(dummy,"%4.1f",farea_max);
            TText *tparameter12=new TText(.7,.15,dummy);
            tparameter12->Draw();
        Tl.DrawLatex(0.05,0.057,"analyzed data record: ");
            sprintf (datafile, "%s%s", filename, ".ASF");
            TText *tparameter13=new TText(.5,.05,datafile);
            tparameter13->Draw();
    coverview->cd(5);//fith window in canvas
        hpos_cut->GetXaxis()->SetTitle("Position x [µm]");
        hpos_cut->GetYaxis()->SetTitle("Position y [µm]");
        gStyle->SetPalette(1);
        hpos_cut->Draw("colz");
    coverview->cd(6);//sixth window in canvas
        coverview_6->SetLogy(1);
        hecc_cut->GetXaxis()->SetTitle("e = b/a [1]");
        hecc_cut->SetLineColor(2);
        hecc_cut->SetLineWidth(2);
        hecc_cut->Draw();
        hecc->SetLineColor(1);
        hecc->SetLineStyle(3);
        hecc->Draw("same");
    coverview->cd(7);//seventh window in canvas
        harea_ca_cut->GetXaxis()->SetTitle("calculated area [µm^{2}]");
        harea_ca_cut->SetLineColor(2);
        harea_ca_cut->SetLineWidth(2);
        harea_ca_cut->Draw();
        fgaus->SetRange(farea_min,farea_max);
        harea_ca_cut->Fit("fgaus","R");
        fgaus->SetLineColor(3);
        harea_ca->SetLineColor(1);
        harea_ca->SetLineStyle(3);
        harea_ca->Draw("same");
        fgaus->Draw("same");
    coverview->cd(8);//eighth window in canvas
        TLatex T2;
        T2.SetTextAlign(12);
        T2.SetTextSize(0.04);
        T2.DrawLatex(.05,.95,"Filterlevels  #rightarrow  Number of filtered craters");
        T2.DrawLatex(.05,.85,"with different cuts without combinations:");
        T2.DrawLatex(0.05,0.757,"1)   Position:");
            sprintf(dummy,"%4.0f",anzahl-3-PositionFilter);
            TText *tparameter14=new TText(.5,.75,dummy);
            tparameter14->Draw();
        T2.DrawLatex(0.05,0.657,"2)   Axis:");
            sprintf(dummy,"%4.0f",anzahl-3-AxisFilter);
            TText *tparameter15=new TText(.5,.65,dummy);
            tparameter15->Draw();
        T2.DrawLatex(0.05,0.557,"3)   Eccentricity:");
            sprintf(dummy,"%4.0f",anzahl-3-EccentricityFilter);
            TText *tparameter16=new TText(.5,.55,dummy);
            tparameter16->Draw();
        T2.DrawLatex(0.05,0.457,"4)   Area:");
            sprintf(dummy,"%4.0f",anzahl-3-AreaFilter);
            TText *tparameter17=new TText(.5,.45,dummy);
            tparameter17->Draw();
        T2.DrawLatex(0.05,0.357,"5)   Underground:");
            sprintf(dummy,"%4.0f",K_ua);
            TText *tparameter18=new TText(.5,.35,dummy);
            tparameter18->Draw();
        T2.DrawLatex(0.05,0.257,"6)   Total:");
            sprintf(dummy,"%4.0f",anzahl-3-Ncrater_filtered);
            TText *tparameter19=new TText(.5,.25,dummy);
            tparameter19->Draw();
        T2.DrawLatex(0.05,0.157,"r_{EA},  r_{EB},  r_{AB}:");
            sprintf(dummy,"%4.3f",r_EA);
            TText *tparameter20=new TText(.3,.15,dummy);
            tparameter20->Draw();
            sprintf(dummy,"%4.3f",r_EB);
            TText *tparameter21=new TText(.5,.15,dummy);
            tparameter21->Draw();
            sprintf(dummy,"%4.3f",r_AB);
            TText *tparameter22=new TText(.7,.15,dummy);
            tparameter22->Draw();
        T2.DrawLatex(0.05,0.057,"A_{Scan} [mm^{2}] = ");
            sprintf(dummy,"%4.3f",A_scan/1000000);
            TText *tparameter23=new TText(.3,.05,dummy);
            tparameter23->Draw();

   //Output - Number Of Filtered Craters:
    cout<< "Number Of Filtered Craters: " << Ncrater_filtered <<endl;

//OUTPUT//END////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

//Saving Histograms Into Files__FUNCTION_DO_SAVE_________________________________________________________________________________________________________________________________
void TCraterAnalysis::DoSave()
{
   //Defining Programtype:
    if (!(fPosCheck->IsOn())){
        sprintf(program,"CraterAnalysis");
    }
    if (fPosCheck->IsOn() && fPosCheckMethod[0]->IsOn()){
        sprintf(program,"CraterAnalysisPosCheckEdge");
    }
    if (fPosCheck->IsOn() && fPosCheckMethod[1]->IsOn()){
        sprintf(program,"CraterAnalysisPosCheckTotal");
    }

   //Defining Path:
    sprintf (datafile, "%s%s%s", path, filename, ".ASF");//path of analyzed data record
    char savefile[200];//path for saved histograms

    ofstream outfile;

   //h_b_axis:
    sprintf (savefile, "%s%s%s%s", path, "h_b_axis_", filename, ".txt");
    outfile.open(savefile);
    //Header:
    outfile<<"Histogram of Semi Minor Axis without Cuts:  h_b_axis"<<"\n"
           <<"analyzed data file:  "<<datafile<<"\n"
           <<"Analysis-Program:  "<<program<<"\n"
           <<"scanned area [Fields x,y]:  "<<image_count_x+1<<"\t"<<image_count_y+1<<"\t"<<"A_Scan [µm²] = "<<"\t"<<A_scan<<"\n"
           <<"Cuts [x_min, x_max, y_min, y_max, b_min, b_max, e_min, e_max, area_min, area_max]:  "<<"\n"
           <<fx_min<<"  "<<fx_max<<"  "<<fy_min<<"  "<<fy_max<<"  "<<fb_min<<"  "<<fb_max<<"  "<<fe_min<<"  "<<fe_max<<"  "<<farea_min<<"  "<<farea_max<<"\n"
           <<"Correlation Coefficients [r_EA, r_EB, r_AB]:  "<<r_EA<<"\t"<<r_EB<<"\t"<<r_AB<<"\n"
           <<"Number of unfiltered Craters = "<<Ncrater_unfiltered<<"\t"<<"Number of filtered Craters = "<<Ncrater_filtered<<"\n"
           <<"Scaling-Factors [x,y]: "<<scalex<<"\t"<<scaley<<"\n"
           <<"Results [Bin Number, Bin Center, Bin Content]:"<<"\n"
           <<"first line: Underflows, last line: Overflows"<<endl;
    //Results:
    for (Int_t i=0;i<h_b_axis->GetNbinsX()+2;i++) {
        outfile<<i<<"\t"<<h_b_axis->GetBinCenter(i)<<"\t"<<h_b_axis->GetBinContent(i)<<endl;
    }
    outfile.close();

   //h_b_axis_cut:
    sprintf (savefile, "%s%s%s%s", path, "h_b_axis_cut_", filename, ".txt");
    outfile.open(savefile);
    //Header:
    outfile<<"Histogram of Semi Minor Axis with Cuts:  h_b_axis_cut"<<"\n"
           <<"analyzed data file:   "<<datafile<<"\n"
           <<"Analysis-Program:  "<<program<<"\n"
           <<"scanned area [Fields x,y]:  "<<image_count_x+1<<"\t"<<image_count_y+1<<"\t"<<"A_Scan [µm²] = "<<"\t"<<A_scan<<"\n"
           <<"Cuts [x_min, x_max, y_min, y_max, b_min, b_max, e_min, e_max, area_min, area_max]:  "<<"\n"
           <<fx_min<<"  "<<fx_max<<"  "<<fy_min<<"  "<<fy_max<<"  "<<fb_min<<"  "<<fb_max<<"  "<<fe_min<<"  "<<fe_max<<"  "<<farea_min<<"  "<<farea_max<<"\n"
           <<"Correlation Coefficients [r_EA, r_EB, r_AB]:  "<<r_EA<<"\t"<<r_EB<<"\t"<<r_AB<<"\n"
           <<"Number of unfiltered Craters = "<<Ncrater_unfiltered<<"\t"<<"Number of filtered Craters = "<<Ncrater_filtered<<"\n"
           <<"Scaling-Factors [x,y]: "<<scalex<<"\t"<<scaley<<"\n"
           <<"Results [Bin Number, Bin Center, Bin Content]:"<<"\n"
           <<"first line: Underflows, last line: Overflows"<<endl;
    //Results:
    for (Int_t i=0;i<h_b_axis_cut->GetNbinsX()+2;i++) {
        outfile<<i<<"\t"<<h_b_axis_cut->GetBinCenter(i)<<"\t"<<h_b_axis_cut->GetBinContent(i)<<endl;
    }
    outfile.close();

   //hecc:
    sprintf (savefile, "%s%s%s%s", path, "hecc_", filename, ".txt");
    outfile.open(savefile);
    //Header:
    outfile<<"Histogram of Eccentricity without Cuts:  hecc"<<"\n"
           <<"analyzed data file:  "<<datafile<<"\n"
           <<"Analysis-Program:  "<<program<<"\n"
           <<"scanned area [Fields x,y]:  "<<image_count_x+1<<"\t"<<image_count_y+1<<"\t"<<"A_Scan [µm²] = "<<"\t"<<A_scan<<"\n"
           <<"Cuts [x_min, x_max, y_min, y_max, b_min, b_max, e_min, e_max, area_min, area_max]:  "<<"\n"
           <<fx_min<<"  "<<fx_max<<"  "<<fy_min<<"  "<<fy_max<<"  "<<fb_min<<"  "<<fb_max<<"  "<<fe_min<<"  "<<fe_max<<"  "<<farea_min<<"  "<<farea_max<<"\n"
           <<"Correlation Coefficients [r_EA, r_EB, r_AB]:  "<<r_EA<<"\t"<<r_EB<<"\t"<<r_AB<<"\n"
           <<"Number of unfiltered Craters = "<<Ncrater_unfiltered<<"\t"<<"Number of filtered Craters = "<<Ncrater_filtered<<"\n"
           <<"Scaling-Factors [x,y]: "<<scalex<<"\t"<<scaley<<"\n"
           <<"Results [Bin Number, Bin Center, Bin Content]:"<<"\n"
           <<"first line: Underflows, last line: Overflows"<<endl;
    //Results:
    for (Int_t i=0;i<hecc->GetNbinsX()+2;i++) {
        outfile<<i<<"\t"<<hecc->GetBinCenter(i)<<"\t"<<hecc->GetBinContent(i)<<endl;
    }
    outfile.close();

   //hecc_cut:
    sprintf (savefile, "%s%s%s%s", path, "hecc_cut_", filename, ".txt");
    outfile.open(savefile);
    //Header:
    outfile<<"Histogram of Eccentricity with Cuts:  hecc_cut"<<"\n"
           <<"analyzed data file:  "<<datafile<<"\n"
           <<"Analysis-Program:  "<<program<<"\n"
           <<"scanned area [Fields x,y]:  "<<image_count_x+1<<"\t"<<image_count_y+1<<"\t"<<"A_Scan [µm²] = "<<"\t"<<A_scan<<"\n"
           <<"Cuts [x_min, x_max, y_min, y_max, b_min, b_max, e_min, e_max, area_min, area_max]:  "<<"\n"
           <<fx_min<<"  "<<fx_max<<"  "<<fy_min<<"  "<<fy_max<<"  "<<fb_min<<"  "<<fb_max<<"  "<<fe_min<<"  "<<fe_max<<"  "<<farea_min<<"  "<<farea_max<<"\n"
           <<"Correlation Coefficients [r_EA, r_EB, r_AB]:  "<<r_EA<<"\t"<<r_EB<<"\t"<<r_AB<<"\n"
           <<"Number of unfiltered Craters = "<<Ncrater_unfiltered<<"\t"<<"Number of filtered Craters = "<<Ncrater_filtered<<"\n"
           <<"Scaling-Factors [x,y]: "<<scalex<<"\t"<<scaley<<"\n"
           <<"Results [Bin Number, Bin Center, Bin Content]:"<<"\n"
           <<"first line: Underflows, last line: Overflows"<<endl;
    //Results:
    for (Int_t i=0;i<hecc_cut->GetNbinsX()+2;i++) {
        outfile<<i<<"\t"<<hecc_cut->GetBinCenter(i)<<"\t"<<hecc_cut->GetBinContent(i)<<endl;
    }
    outfile.close();

   //harea_ea:
    sprintf (savefile, "%s%s%s%s", path, "harea_ea_", filename, ".txt");
    outfile.open(savefile);
    //Header:
    outfile<<"Histogram of Enclosed Area without Cuts:  harea_ea"<<"\n"
           <<"analyzed data file:  "<<datafile<<"\n"
           <<"Analysis-Program:  "<<program<<"\n"
           <<"scanned area [Fields x,y]:  "<<image_count_x+1<<"\t"<<image_count_y+1<<"\t"<<"A_Scan [µm²] = "<<"\t"<<A_scan<<"\n"
           <<"Cuts [x_min, x_max, y_min, y_max, b_min, b_max, e_min, e_max, area_min, area_max]:  "<<"\n"
           <<fx_min<<"  "<<fx_max<<"  "<<fy_min<<"  "<<fy_max<<"  "<<fb_min<<"  "<<fb_max<<"  "<<fe_min<<"  "<<fe_max<<"  "<<farea_min<<"  "<<farea_max<<"\n"
           <<"Correlation Coefficients [r_EA, r_EB, r_AB]:  "<<r_EA<<"\t"<<r_EB<<"\t"<<r_AB<<"\n"
           <<"Number of unfiltered Craters = "<<Ncrater_unfiltered<<"\t"<<"Number of filtered Craters = "<<Ncrater_filtered<<"\n"
           <<"Scaling-Factors [x,y]: "<<scalex<<"\t"<<scaley<<"\n"
           <<"Results [Bin Number, Bin Center, Bin Content]:"<<"\n"
           <<"first line: Underflows, last line: Overflows"<<endl;
    //Results:
    for (Int_t i=0;i<harea_ea->GetNbinsX()+2;i++) {
        outfile<<i<<"\t"<<harea_ea->GetBinCenter(i)<<"\t"<<harea_ea->GetBinContent(i)<<endl;
    }
    outfile.close();

   //harea_ca:
    sprintf (savefile, "%s%s%s%s", path, "harea_ca_", filename, ".txt");
    outfile.open(savefile);
    //Header:
    outfile<<"Histogram of Calculated Area without Cuts:  harea_ca"<<"\n"
           <<"analyzed data file:  "<<datafile<<"\n"
           <<"Analysis-Program:  "<<program<<"\n"
           <<"scanned area [Fields x,y]:  "<<image_count_x+1<<"\t"<<image_count_y+1<<"\t"<<"A_Scan [µm²] = "<<"\t"<<A_scan<<"\n"
           <<"Cuts [x_min, x_max, y_min, y_max, b_min, b_max, e_min, e_max, area_min, area_max]:  "<<"\n"
           <<fx_min<<"  "<<fx_max<<"  "<<fy_min<<"  "<<fy_max<<"  "<<fb_min<<"  "<<fb_max<<"  "<<fe_min<<"  "<<fe_max<<"  "<<farea_min<<"  "<<farea_max<<"\n"
           <<"Correlation Coefficients [r_EA, r_EB, r_AB]:  "<<r_EA<<"\t"<<r_EB<<"\t"<<r_AB<<"\n"
           <<"Number of unfiltered Craters = "<<Ncrater_unfiltered<<"\t"<<"Number of filtered Craters = "<<Ncrater_filtered<<"\n"
           <<"Scaling-Factors [x,y]: "<<scalex<<"\t"<<scaley<<"\n"
           <<"Results [Bin Number, Bin Center, Bin Content]:"<<"\n"
           <<"first line: Underflows, last line: Overflows"<<endl;
    //Results:
    for (Int_t i=0;i<harea_ca->GetNbinsX()+2;i++) {
        outfile<<i<<"\t"<<harea_ca->GetBinCenter(i)<<"\t"<<harea_ca->GetBinContent(i)<<endl;
    }
    outfile.close();

   //harea_ca_cut:
    sprintf (savefile, "%s%s%s%s", path, "harea_ca_cut_", filename, ".txt");
    outfile.open(savefile);
    //Header:
    outfile<<"Histogram of Calculated Area with Cuts:  harea_ca_cut"<<"\n"
           <<"analyzed data file:   "<<datafile<<"\n"
           <<"Analysis-Program:  "<<program<<"\n"
           <<"scanned area [Fields x,y]:  "<<image_count_x+1<<"\t"<<image_count_y+1<<"\t"<<"A_Scan [µm²] = "<<"\t"<<A_scan<<"\n"
           <<"Cuts [x_min, x_max, y_min, y_max, b_min, b_max, e_min, e_max, area_min, area_max]:  "<<"\n"
           <<fx_min<<"  "<<fx_max<<"  "<<fy_min<<"  "<<fy_max<<"  "<<fb_min<<"  "<<fb_max<<"  "<<fe_min<<"  "<<fe_max<<"  "<<farea_min<<"  "<<farea_max<<"\n"
           <<"Correlation Coefficients [r_EA, r_EB, r_AB]:  "<<r_EA<<"\t"<<r_EB<<"\t"<<r_AB<<"\n"
           <<"Number of unfiltered Craters = "<<Ncrater_unfiltered<<"\t"<<"Number of filtered Craters = "<<Ncrater_filtered<<"\n"
           <<"Scaling-Factors [x,y]: "<<scalex<<"\t"<<scaley<<"\n"
           <<"Results [Bin Number, Bin Center, Bin Content]:"<<"\n"
           <<"first line: Underflows, last line: Overflows"<<endl;
    //Results:
    for (Int_t i=0;i<harea_ca_cut->GetNbinsX()+2;i++) {
        outfile<<i<<"\t"<<harea_ca_cut->GetBinCenter(i)<<"\t"<<harea_ca_cut->GetBinContent(i)<<endl;
    }
    outfile.close();

    cout<<"DoSave"<<endl;
}

//Retrieving The Main Funtion__FUNCTION_CRATER_ANALYSIS__________________________________________________________________________________________________________________________
void CraterAnalysis()
{
    new TCraterAnalysis();
}
