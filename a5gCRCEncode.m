function [all_data] = a5gCRCEncode(input, type)
%% inputlar:
% input: CRC bitleri eklenecek olan input bitleri
% type: CRC tipini belirtir (string - '24A'-'24B'-'24C'-'16'-'11'-'6')

%% output: all_data
if strcmpi(type, '24A')
    gx = zeros(1,25);
    %mert
    index = 25-[24 23 18 17 14 11 10 7 6 5 4 3 1 0];
    % index = [25 24 19 18 15 12 11 8 7 6 5 4 2 1];
    gx(index) = 1;
elseif strcmpi(type, '24B')
    gx = zeros(1,25);
    %mert
    index = 25-[24 23 6 5 1 0];
    % index = [25 24 7 6 2 1];
    gx(index) = 1;
elseif strcmpi(type, '24C')
    gx = zeros(1,25);
    %mert
    index = 25-[24 23 21 20 17 15 13 12 8 4 2 1 0];
    % index = [25 24 22 21 18 16 14 13 9 5 3 2 1];
    gx(index) = 1;
elseif strcmpi(type, '16')
    gx = zeros(1,17);
    %mert
    index = 17-[16 12 5 0];
    % index = [17 13 6 1];
    gx(index) = 1;
elseif strcmpi(type, '13')
    gx = zeros(1,14);
    index = 14- [13 12 10  8 6 5 4 3 2 1 0];
   % index = 13-[11 10 9 8 4 1 0];
    gx(index) = 1;       
elseif strcmpi(type, '12')
    gx = zeros(1,13);
    index = 13- [12 9 8 3 2 1 0];
   % index = 13-[11 10 9 8 4 1 0];
    %index = 13-[12 9 5 4 3 1 0];%123B
    % index = 13-[12 10 7 5 4 3 2]; %1957
   % index = 13 -[12 10  7 4 3 2 1]; %fasura
    gx(index) = 1;    
elseif strcmpi(type, '11')
    gx = zeros(1,12);
    %mert
    index = 12-[11 10 9 5 0];
    %index = 12-[11 10 9 3 2 1 0]; %a0f
    %index = 12-[11 8 5 4 3 2 1 0]; %93f
   % index = 12-[11 9 7 4 3 2 0]; %a9d
    gx(index) = 1;
elseif strcmpi(type, '6')
    gx = zeros(1,7);
    %mert
    index = 7-[6 5 0];
    % index = [7 6 1];
    gx(index) = 1;
else
    %warning('Invalid CRC type or input length!');
    disp('Invalid CRC type or input length!');
end

px = input;
pol_length=length(gx)-1;
pxr= [px zeros(1,pol_length)];

yy = [0, pxr(1:pol_length)];
for i=1:length(px)
    a = [yy(2:end), pxr(i+pol_length)];
    if a(1) == 1
        yy = rem(gx+a,2);
    else
        yy = a;
    end
end
all_data = [input yy(2:end)];


end