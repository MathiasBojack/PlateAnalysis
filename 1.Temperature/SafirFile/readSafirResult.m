function readSafirResult(h)

fReadName  = ['Heat_Transfer_' num2str(h*100) 'cm.out'];
fRid       = fopen(fReadName,'r');

fWriteName = [num2str(h*100) 'cm.txt'];
fWid       = fopen(fWriteName,'W');
count = 0;
while ~feof(fRid)

    tline = fgetl(fRid);
    if contains(tline,'TIME = ')
        count = count + 1;
        if count == 601
            break
        end
        time   = regexp(tline,'\d*','Match');  % cell type
        tline  = fgetl(fRid);
        tline  = fgetl(fRid);
        tline  = fgetl(fRid);
        tline  = fgetl(fRid);
        tline  = fgetl(fRid);
        fRead  = fscanf(fRid,'%d %f', [2, 46]);
        fprintf(fWid,'%8.0f',str2num(time{1}));
        for i = 1:length(fRead(2,:))
         fprintf(fWid,'%8.1f',fRead(2,i));
        end
        fprintf(fWid,'%s\n', '');
        
    end
    
end
fclose(fRid);
fclose(fWid);
end