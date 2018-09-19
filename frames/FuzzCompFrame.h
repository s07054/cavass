/*
  Copyright 1993-2013, 2015 Medical Image Processing Group
              Department of Radiology
            University of Pennsylvania

This file is part of CAVASS.

CAVASS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAVASS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CAVASS.  If not, see <http://www.gnu.org/licenses/>.

*/

//======================================================================
/**
 * \file   FuzzCompFrame.h
 * \brief  FuzzCompFrame definition.
 * \author Andre Souza, Ph.D.
 *
 * Copyright: (C) 2007, CAVASS
 *
 * Rise and shine and give God your glory (glory).
 */
//======================================================================
#ifndef __FuzzCompFrame_h
#define __FuzzCompFrame_h

#include  <algorithm>
#include  <stack>
#include  "wx/dnd.h"
#include  "wx/docview.h"
#include  "wx/splitter.h"
#include  "MainFrame.h"

#ifndef WIN32
    #include  <unistd.h>
#endif
#include  <stdlib.h>

class  CineControls;
class  ColorControls;
class  GrayMapControls;
class  SaveScreenControls;
class  SetFuzzCompIndexControls;
class  ClassifyControls;
class  HistoSettingControls;
class  CovarianceControls;
class  ParameterControls;
class  FunctionControls;
class  FuzzCompCanvas;

/** \brief FuzzCompFrame class definition.
 *
 *  a frame for the fuzzy connectedness segmentation method.
 */
class FuzzCompFrame : public MainFrame {
    CineControls*              mCineControls;
    GrayMapControls*           mGrayMap1Controls;
    GrayMapControls*           mGrayMap2Controls;
    GrayMapControls*           mGrayMap3Controls;
    GrayMapControls*           mGrayMap4Controls;
    SaveScreenControls*        mSaveScreenControls;
    SetFuzzCompIndexControls*  mSetIndex1Controls;    
    ClassifyControls*          mClassifyControls;
    HistoSettingControls*      mHistoSettingControls;
    CovarianceControls*        mCovarianceControls;
    ParameterControls*         mParameterControls;
    FunctionControls*          mFunctionControls;

  public:
    static wxFileHistory  _fileHistory;
	FuzzCompCanvas*        mFuzzCanvas;  
    enum {
        ID_PREVIOUS=ID_LAST, ID_NEXT,
        ID_SET_INDEX1, ID_MATCH_INDEX1,
        ID_GRAYMAP1,         ID_GRAYMAP2,
        ID_GRAYMAP3,         ID_GRAYMAP4,        
        ID_CT_LUNG,        ID_CT_SOFT_TISSUE, ID_CT_BONE, ID_PET,
        ID_APPLYTRAIN, ID_CLEARSEED,
        ID_CINE,             ID_RESET,
        ID_CINE_TIMER, ID_CINE_SLIDER, ID_CINE_FORWARD, ID_CINE_FORWARD_BACKWARD,
        ID_CENTER1_SLIDER, ID_CENTER2_SLIDER, ID_CENTER3_SLIDER, ID_CENTER4_SLIDER,
        ID_WIDTH1_SLIDER,  ID_WIDTH2_SLIDER, ID_WIDTH3_SLIDER, ID_WIDTH4_SLIDER,
        ID_INVERT1,        ID_INVERT2, ID_INVERT3, ID_INVERT4,
        ID_SLICE1_SLIDER,
        ID_SCALE1_SLIDER,
        ID_ZOOM_IN, ID_ZOOM_OUT,
        ID_OVERLAY, ID_BUSY_TIMER,
        ID_OVERWRITE_SCREEN, ID_APPEND_SCREEN, ID_BROWSE_SCREEN,
        ID_PROCESS, ID_SAVE
    };
    enum {
        WHICH_BOTH, WHICH_FIRST, WHICH_SECOND
    };
    enum {
        ID_AFFIN_COMBOBOX, ID_MODE_COMBOBOX, ID_ADJ_COMBOBOX,ID_BRUSH_SLIDER, 
        ID_THRES_SLIDER,ID_PARALLEL_INDEX, ID_LOWBIN_SLIDER, ID_HIGHBIN_SLIDER,
        ID_HIGH, ID_LOW, ID_DIFF, ID_SUM, ID_RELDIFF, ID_HIGH_COMBO, ID_LOW_COMBO,
        ID_DIFF_COMBO, ID_SUM_COMBO, ID_RELDIFF_COMBO, ID_AFFOUT_COMBO,
		ID_STARTSLICE_SLIDER, ID_ENDSLICE_SLIDER,
        ID_WEIGHT1_SLIDER, ID_LEVEL1_SLIDER,ID_WIDTHLEVEL1_SLIDER,
        ID_WEIGHT2_SLIDER, ID_LEVEL2_SLIDER,ID_WIDTHLEVEL2_SLIDER,
        ID_WEIGHT3_SLIDER, ID_LEVEL3_SLIDER,ID_WIDTHLEVEL3_SLIDER,
        ID_WEIGHT4_SLIDER, ID_LEVEL4_SLIDER,ID_WIDTHLEVEL4_SLIDER,
        ID_WEIGHT5_SLIDER, ID_LEVEL5_SLIDER,ID_WIDTHLEVEL5_SLIDER
    };
    int            mFileOrDataCount;    ///< two needed for overlay
    
//    wxSlider*      m_center;            ///< contrast center slider
//    wxSlider*      m_width;             ///< contrast width slider
  protected:
    //wxString       mTargetName;
    wxString       mSourceName, mTmpName;
    wxButton*      mSave;
	wxButton*      mApplyTrain;
    wxButton*      mClearSeed;
	wxButton*      mProcess;
    wxMenu*        m_options_menu;      ///< options menu
    wxSizer*       mBottomSizer;
    wxStaticBox*   mSetIndex1Box;
    wxSizer*       mSetIndex1Sizer;
    wxStaticBox*   mClassifyBox;
    wxSizer*       mClassifySizer;

    //cine related items
    bool           mForward, mForwardBackward, mDirectionIsForward;  ///< cine direction
    wxTimer*       m_cine_timer;        ///< time for cine

    void initializeMenu ( void );
    void addButtonBox ( void );

    void OnCineTimer           ( wxTimerEvent& e    );
    void OnCineForward         ( wxCommandEvent& e  );
    void OnCineForwardBackward ( wxCommandEvent& e  );
    void OnCineSlider          ( wxScrollEvent& e   );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  public:
    static void createFuzzCompFrame ( wxFrame* parentFrame, bool useHistory=true );
    FuzzCompFrame ( bool maximize=false, int w=800, int h=600 );
    ~FuzzCompFrame ( void );
	static bool match ( wxString filename );
	virtual bool filenameMatch ( wxString filename ) const {
	    return match( filename );
	};

	virtual void OnInput ( wxCommandEvent& unused );
    void loadFile ( const char* const fname );
    void loadData ( char* name,
        const int xSize, const int ySize, const int zSize,
        const double xSpacing, const double ySpacing, const double zSpacing,
        const int* const data,
        const ViewnixHeader* const vh=NULL, const bool vh_initialized=false );

    void OnCenter1Slider ( wxScrollEvent& e );    
    void OnWidth1Slider  ( wxScrollEvent& e );    
    void OnSlice1Slider  ( wxScrollEvent& e );    
    void OnScale1Slider  ( wxScrollEvent& e );

    void OnCenter2Slider ( wxScrollEvent& e );    
    void OnWidth2Slider  ( wxScrollEvent& e );    
    
    
    void OnCenter3Slider ( wxScrollEvent& e );    
    void OnWidth3Slider  ( wxScrollEvent& e );    
    

    void OnCenter4Slider ( wxScrollEvent& e );    
    void OnWidth4Slider  ( wxScrollEvent& e );    

	void OnCTLung       ( wxCommandEvent& unused );
	void OnCTSoftTissue ( wxCommandEvent& unused );
	void OnCTBone       ( wxCommandEvent& unused );
    void OnPET          ( wxCommandEvent& unused );

    void OnBusyTimer     ( wxTimerEvent& e  );
#ifdef __WXX11__
    void OnUpdateUICenter1Slider ( wxUpdateUIEvent& e );
    void OnUpdateUIWidth1Slider  ( wxUpdateUIEvent& e );
    void OnUpdateUISlice1Slider  ( wxUpdateUIEvent& e );
    void OnUpdateUIScale1Slider  ( wxUpdateUIEvent& e );

    void OnUpdateUICenter2Slider ( wxUpdateUIEvent& e );
    void OnUpdateUIWidth2Slider  ( wxUpdateUIEvent& e );
    
    void OnUpdateUICenter3Slider ( wxUpdateUIEvent& e );
    void OnUpdateUIWidth3Slider  ( wxUpdateUIEvent& e );
    
    void OnUpdateUICenter4Slider ( wxUpdateUIEvent& e );
    void OnUpdateUIWidth4Slider  ( wxUpdateUIEvent& e );
    

    void OnUpdateUICineSlider    ( wxUpdateUIEvent& e );
#endif
    void OnPrevious     ( wxCommandEvent& unused );
    void OnNext         ( wxCommandEvent& unused );
    void OnSetIndex1    ( wxCommandEvent& unused );
    void OnGrayMap1     ( wxCommandEvent& unused );
    void OnGrayMap2     ( wxCommandEvent& unused );
    void OnGrayMap3     ( wxCommandEvent& unused );
    void OnGrayMap4     ( wxCommandEvent& unused );
    void OnCine         ( wxCommandEvent& unused );
    void OnReset        ( wxCommandEvent& unused );
    void OnProcess      ( wxCommandEvent& unused );
    void OnSave         ( wxCommandEvent& unused );
    void OnMode         ( wxCommandEvent& unused );   
	void OnApplyTrain    ( wxCommandEvent& unused );               
    void OnClearSeed    ( wxCommandEvent& unused );               
    void OnThresSlider  ( wxScrollEvent& e );
    void OnBrushSlider  ( wxScrollEvent& e );
    void OnLowBinSlider  ( wxScrollEvent& e );
    void OnHighBinSlider  ( wxScrollEvent& e );
    void OnAffinity     ( wxCommandEvent& unused );
	void OnHigh         ( wxCommandEvent& e );
	void OnLow         ( wxCommandEvent& e );
	void OnDiff         ( wxCommandEvent& e );
	void OnSum         ( wxCommandEvent& e );
	void OnRelDiff         ( wxCommandEvent& e );
    void OnHighCombo    ( wxCommandEvent& unused );
    void OnLowCombo    ( wxCommandEvent& unused );
    void OnDiffCombo    ( wxCommandEvent& unused );
    void OnSumCombo    ( wxCommandEvent& unused );
    void OnRelDiffCombo    ( wxCommandEvent& unused );
	void OnAdjacency   ( wxCommandEvent& e );
	void OnPararell    ( wxCommandEvent& e );
    
	void OnSeed    ( void );
	void OffSave   ( void );

    void OnWeight1Slider ( wxScrollEvent& e );
    void OnLevel1Slider ( wxScrollEvent& e );
    void OnWidthLevel1Slider ( wxScrollEvent& e );

    void OnWeight2Slider ( wxScrollEvent& e );
    void OnLevel2Slider ( wxScrollEvent& e );
    void OnWidthLevel2Slider ( wxScrollEvent& e );

    void OnWeight3Slider ( wxScrollEvent& e );
    void OnLevel3Slider ( wxScrollEvent& e );
    void OnWidthLevel3Slider ( wxScrollEvent& e );

    void OnWeight4Slider ( wxScrollEvent& e );
    void OnLevel4Slider ( wxScrollEvent& e );
    void OnWidthLevel4Slider ( wxScrollEvent& e );

    void OnWeight5Slider ( wxScrollEvent& e );
    void OnLevel5Slider ( wxScrollEvent& e );
    void OnWidthLevel5Slider ( wxScrollEvent& e );
	
	void OnStartSliceSlider ( wxScrollEvent& e );
	void OnEndSliceSlider ( wxScrollEvent& e );


    void OnOverwriteScreen ( wxCommandEvent& unused );
    void OnAppendScreen    ( wxCommandEvent& unused );
    void OnBrowseScreen    ( wxCommandEvent& unused );

    void OnSaveScreen   ( wxCommandEvent& e );
    void OnCopy         ( wxCommandEvent& e );
    void OnOpen         ( wxCommandEvent& e );
    void OnZoomIn       ( wxCommandEvent& e );
    void OnZoomOut      ( wxCommandEvent& e );
    void OnOverlay      ( wxCommandEvent& e );
    void OnInvert1      ( wxCommandEvent& e );
    void OnInvert2      ( wxCommandEvent& e );
    void OnInvert3      ( wxCommandEvent& e );
    void OnInvert4      ( wxCommandEvent& e );
    void OnHideControls ( wxCommandEvent& e );
    void OnPrintPreview ( wxCommandEvent& e );
    void OnPrint        ( wxCommandEvent& e );
    void OnChar         ( wxKeyEvent& e );
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    DECLARE_DYNAMIC_CLASS( FuzzCompFrame )
    DECLARE_EVENT_TABLE()
};

#endif
//======================================================================
