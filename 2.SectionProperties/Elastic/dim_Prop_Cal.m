function Input = dim_Prop_Cal( Input)

Adim_Prop = adim_Prop_Cal( Input.Load.Time,Input.Geometry.h,Input.Properties.Nu );

A0 = Adim_Prop.A0;
B0 = Adim_Prop.B0;
D0 = Adim_Prop.D0;
NT0 = Adim_Prop.NT0;
MT0 = Adim_Prop.MT0;

Input.Properties.De = Input.Properties.E/(1-Input.Properties.Nu^2)*Input.Geometry.h^3*(A0*D0-B0^2)/A0;
Input.Properties.CHI_SS = Adim_Prop.CHI_SS;   % thermal curvature for simply supported boundary
Input.Properties.CHI_FB = Adim_Prop.CHI_FB;   % thermal curvature for free boundary
Input.Properties.A = Input.Properties.E*Input.Geometry.h/(1-Input.Properties.Nu^2)*A0;
Input.Properties.B = Input.Properties.E*Input.Geometry.h^2/(1-Input.Properties.Nu^2)*B0;
Input.Properties.D = Input.Properties.E*Input.Geometry.h^3/(1-Input.Properties.Nu^2)*D0;
Input.Properties.NT = -Input.Properties.E*Input.Geometry.h/(1-Input.Properties.Nu)*NT0;
Input.Properties.MT = Input.Properties.E*Input.Geometry.h^2/(1-Input.Properties.Nu)*MT0;