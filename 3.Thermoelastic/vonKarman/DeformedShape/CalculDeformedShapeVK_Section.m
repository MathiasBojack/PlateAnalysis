clear
Height    = 12;
lambda    = 4;
% Thickness = [0.1 0.15 0.2 0.25 0.3];
Thickness = 0.15;
TimeList  = [60 1200:1200:14400];




for i = 1:length(Height)
    for j = 1:length(lambda)
        for k = 1:length(Thickness)
            for l = 1:length(TimeList)
                cnt = (i-1)*length(lambda)*length(Thickness)*length(TimeList)...
                    + (j-1)*length(Thickness)*length(TimeList) ...
                    + (k-1)*length(TimeList) +l;
                variable(:,cnt) = [Height(i) Height(i)/lambda(j) ...
                    Thickness(k) TimeList(l) ]';
            end
        end
    end
end

%%

CNT = length(Height)*length(lambda)*length(Thickness)*length(TimeList);
for cnt = 1:CNT
    Input.Geometry.a    = variable(1,cnt);
    Input.Geometry.b    = variable(2,cnt);
    Input.Geometry.h    = variable(3,cnt);
    Input.Load.Time     = variable(4,cnt);
    Input.Properties.Nu = 0.2;
    Input.Properties.E  = 19.2e9;

    Input.Load.gamma    =  2500*9.8*Input.Geometry.h;
    Input.SolParam.M    = 20;
    Input.SolParam.N    = 20;

    Input = dim_Prop_Cal( Input);

    Solution = Four_ss_rectangular_Plate_VK( Input );
    dname = strcat('./3.Thermoelastic/vonKarman/DeformedShape/Deformee_a_',...
        num2str(Input.Geometry.a),'_b_',num2str(Input.Geometry.b),'_h_',...
        num2str(Input.Geometry.h),'_Time_',num2str(Input.Load.Time),'.mat');
    save(dname,'Solution')
    disp(cnt/CNT)
    
end