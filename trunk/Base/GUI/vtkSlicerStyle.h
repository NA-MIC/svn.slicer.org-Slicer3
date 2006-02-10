#ifndef __vtkSlicerStyle_h
#define __vtkSlicerStyle_h

#include "vtkObject.h"

// Description:
// Definition of Slicer's look and feel.
//
class vtkSlicerStyle : public vtkObject
{
 public:
    static vtkSlicerStyle* New ( );
    vtkTypeRevisionMacro ( vtkSlicerStyle, vtkObject );

    // Description:
    // This method configures the Slicer brand.
    // New style classes can be derived from this one,
    // and the ApplyPresentation method can be overridden
    // to specify different look and feel.
    // Can we plan to adopt Tk tiles to create
    // a platform invariant expression of Slicer brand?
    void ApplyPresentation ( );
    // To be defined later, reading slicer style from xml fikle.
    // for now, just hardcode them.
    void ParseStyleParameters ( );

    // Interface to standard Tk options
    // font and text
    virtual void SetBigFont ( char *font );
    virtual char* GetBigFont ( ) { return this->BigFont ;}
    virtual void SetMedFont ( char *font );
    virtual char* GetMedFont ( ) { return this->MedFont ;}
    virtual void SetSmallFont ( char *font );
    virtual char* GetSmallFont ( ) { return this->SmallFont ;}
    virtual void SetJustify ( char * );
    virtual char* GetJustify ( ) { return this->Justify ;}
    virtual void SetTextLeftJustify ( );
    virtual void SetTextRightJustify ( );
    virtual void SetTextCenterJustify ( );

    // relief
    virtual void SetRelief ( char * );
    virtual void SetFlatRelief ( );
    virtual void SetGrooveRelief ( );

    // color
    // Description:
    // Used by all color-setting methods to do error checking, setting.
    //
    virtual int SetColor ( double *color, double r, double g, double b );
    virtual int SetBgColor ( double *color );
    virtual int SetBgColor ( double r, double g, double b );
    virtual double* GetBgColor ( );
    virtual int SetActiveBgColor ( double r, double g, double b );
    virtual int SetActiveBgColor ( double *color );
    virtual double* GetActiveBgColor ( );
    virtual int SetInsertBgColor ( double r, double g, double b );
    virtual int SetInsertBgColor ( double *color );
    virtual double* GetInsertBgColor ( );
    virtual int SetSelectBgColor ( double r, double g, double b );
    virtual int SetSelectBgColor ( double *color );
    virtual double* GetSelectBgColor ( );
    virtual int SetFgColor ( double r, double g, double b );
    virtual int SetFgColor ( double *color );
    virtual double* GetFgColor ( );
    virtual int SetActiveFgColor ( double r, double g, double b );
    virtual int SetActiveFgColor ( double *color );
    virtual double* GetActiveFgColor ( );
    virtual int SetDisabledFgColor ( double r, double g, double b );
    virtual int SetDisabledFgColor ( double *color );
    virtual double* GetDisabledFgColor ( );
    virtual int SetSelectFgColor ( double r, double g, double b );
    virtual int SetSelectFgColor ( double *color );
    virtual double* GetSelectFgColor ( );
    virtual int SetTroughColor ( double r, double g, double b );
    virtual int SetTroughColor ( double *color );
    virtual double* GetTroughColor ( );

    // highlights
    virtual int SetHighLightColor ( double r, double g, double b );
    virtual int SetHighLightColor ( double *color );
    virtual double* GetHighLightColor ( );
    virtual int SetHighLightBgColor ( double r, double g, double b );
    virtual int SetHighLightBgColor ( double *color );
    virtual double* GetHighLightBgColor ( );
    virtual void SetHighLightThickness ( int thickness );
    virtual int GetHighLightThickness ( ) { return this->HighLightThickness; }
    
    // spacing functions
    virtual void SetBorderWidth ( int width );
    virtual int GetBorderWidth ( ) { return this->BorderWidth ;}
    virtual void SetActiveBorderWidth ( int width );
    virtual int GetActiveBorderWidth ( ) { return this->ActiveBorderWidth ;}
    virtual void SetSelectBorderWidth ( int width );
    virtual int GetSelectBorderWidth ( ) { return this->SelectBorderWidth ;}
    virtual void SetPadX ( int padx );
    virtual int GetPadX ( ) { return this->PadX ;}
    virtual void SetPadY (int pady );
    virtual int GetPadY ( ) { return this->PadY ;}
    
    
 protected:
    // Font and text
    char *BigFont;
    char *MedFont;
    char *SmallFont;
    char *Justify;
    
    // Relief
    char *FlatRelief;
    char *GrooveRelief;
    char *Relief;
    
    // Code colors
    double SavedDataTextColor [3];
    double UnsavedDataTextColor [3];
    double ErrorTextColor [3];
    double ErrorColor [3];
    double WarningTextColor [3];
    double WarningColor [3];
    
    // Foreground and Background color
    double BgColor [3];
    double ActiveBgColor [3];
    double InsertBgColor [3];
    double SelectBgColor [3];
    double FgColor [3];
    double ActiveFgColor [3];
    double DisabledFgColor [3];
    double SelectFgColor [3];
    double TroughColor [3];
    
    // Highlights
    double HighLightColor [3];
    double HighLightBgColor [3];
    int HighLightThickness;
    
    // Spacing
    int BorderWidth;
    int ActiveBorderWidth;
    int SelectBorderWidth;
    int PadX;
    int PadY;

    vtkSlicerStyle ( );
    ~vtkSlicerStyle ( );
        
 private:
    vtkSlicerStyle ( const vtkSlicerStyle&); // Not implemented
    void operator = ( const vtkSlicerStyle&); // Not implemented
};

#endif
