function [ Profil ] = read_Temperature( h )
%READ_TEMPERATURE read temperature distribution  throught the thickness
%   h: thickness of the plate =[10 15 20 25 30] [cm]
%   Profil(1,:) Time
%   Profil(2:end,:) Temperature
fileName = strcat(num2str(h*100),'cm.txt');
fileID = fopen(fileName,'r');
formatSpec = '%f';
sizeT = [47 Inf];

Profil.Temperature = fscanf(fileID,formatSpec,sizeT);
Profil.Position = [-1/2:1/60:0 1/30:1/30:1/2]'*h;
fclose(fileID);
end

