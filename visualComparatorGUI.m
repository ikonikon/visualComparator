function varargout = visualComparatorGUI(varargin)

% Ioannis Konstantelos <ik303@ic.ac.uk>
% Imperial College London
% Version: 0.1
% Created: September 2014
% Modified: November 2014

% VISUALCOMPARATORGUI MATLAB code for visualComparatorGUI.fig
%      VISUALCOMPARATORGUI, by itself, creates a new VISUALCOMPARATORGUI or raises the existing
%      singleton*.
%
%      H = VISUALCOMPARATORGUI returns the handle to a new VISUALCOMPARATORGUI or the handle to
%      the existing singleton*.
%
%      VISUALCOMPARATORGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALCOMPARATORGUI.M with the given input arguments.
%
%      VISUALCOMPARATORGUI('Property','Value',...) creates a new VISUALCOMPARATORGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visualComparatorGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visualComparatorGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visualComparatorGUI

% Last Modified by GUIDE v2.5 18-Nov-2014 12:58:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @visualComparatorGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @visualComparatorGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before visualComparatorGUI is made visible.
function visualComparatorGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visualComparatorGUI (see VARARGIN)

% Choose default command line output for visualComparatorGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes visualComparatorGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = visualComparatorGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in runButton.
function runButton_Callback(hObject, eventdata, handles)
% hObject    handle to runButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Clear Axes
cla(handles.axes1);
cla(handles.axes2);
cla(handles.axes3);
cla(handles.axes4);
cla(handles.axes5);
cla(handles.axes6);

%% Import data and check input validity

try
    Z = importdata(handles.fileName1TextBox.String);
catch
    errordlg('PLease specify valid name for File 1','Error');
end

try
    Q = importdata(handles.fileName2TextBox.String);
catch
    errordlg('PLease specify valid name for File 2','Error');
end

[nObs1,nVar1] = size(Z);
[nObs2,nVar2] = size(Q);

% Display dimensionality
set(handles.nObs1TextBox, 'String', nObs1)
set(handles.nVar1TextBox, 'String', nVar1)
set(handles.nObs2TextBox, 'String', nObs2)
set(handles.nVar2TextBox, 'String', nVar2)

% Error message if nVar1 =~ nVar2
if nVar1 ~= nVar2
    errordlg('Number of  variables does not match','Error');
end

% Determine radio button selection (default/custom)
isDefault = get(handles.radiobuttonDefault, 'Value');

if isDefault == 1
    % define the two sets (half and half)
    set1 = 1:round(nVar1/2);
    set2 = round(nVar1/2) + 1:nVar1;
else
    set1String = get(handles.set1TextBox, 'String');
    set2String = get(handles.set2TextBox, 'String');
    if isempty(set1String) || isempty(set1String)
        % define the two sets (half and half)
        set1 = 1:round(nVar1/2);
        set2 = round(nVar1/2) + 1:nVar1;
    else
        % 
        set1 = eval(set1String);
        set2 = eval(set2String);
    end       
end

% 
if isDefault == 1
    set1String = [num2str(set1(1)),':',num2str(set1(end))];
    set2String = [num2str(set2(1)),':',num2str(set2(end))];
else
    set1String = num2str(set1(1));
    set2String = num2str(set2(1));
end

set(handles.set1StaticBox,'String',set1String)
set(handles.set2StaticBox,'String',set2String)

% Define aggregated variables
Z1 = sum(Z(:,set1),2);
Z2 = sum(Z(:,set2),2);
Q1 = sum(Q(:,set1),2);
Q2 = sum(Q(:,set2),2);

% Check if points are stationary
if all(Z1 == Z1(1))
    errordlg('File 1 : Set 1 is static','Error','modal');
else
    if all(Z2 == Z2(1))
        errordlg('File 1 : Set 2 is static','Error','modal');
    else
        if all(Q1 == Q1(1))
            errordlg('File 2 : Set 1 is static','Error','modal');
        else
            if all(Q2 == Q2(1))
                errordlg('File 2 : Set 2 is static','Error','modal');
            end
        end
    end
end

% Calculate ecdf curve points | f : probability, v : value 
[f_Z1,v_Z1]  = ecdf(Z1);
[f_Z2,v_Z2]  = ecdf(Z2);
[f_Q1,v_Q1]  = ecdf(Q1);
[f_Q2,v_Q2]  = ecdf(Q2);

% define x-axis limits
f_xmin = 0;  
f_xmax = 1;  
v_xmin = min([v_Z1 ; v_Q1]);
v_xmax = max([v_Z1 ; v_Q1]);

% define y-axis limits
f_ymin = 0;
f_ymax = 1;
v_ymin = min([v_Z2 ; v_Q2]);
v_ymax = max([v_Z2 ; v_Q2]);

% adjust axis-limits for stationary variables
if v_xmin == v_xmax
    v_xmin = v_xmin - 1;
    v_xmax = v_xmax + 1;
end
if v_ymin == v_ymax
    v_ymin = v_ymin - 1;
    v_ymax = v_ymax + 1;
end

%% Compute Kolmogorov-Smirnov metrics

% Perform K-S tests
[~,p1,~] = kstest2(Z1,Q1);
[~,p2,~] = kstest2(Z2,Q2);

d1 = nan(size(v_Q1,1),1);
% Find location of maximum ecdf discrepancy for set 1
for i = 1:size(v_Q1,1)
    % find closest point of Z1 to Q1
    [~,idx1] = min(abs(v_Z1-v_Q1(i)));
    d1(i) = f_Z1(idx1)-f_Q1(i);
end

d2 = nan(size(v_Q2,1),1);
% Find location of maximum ecdf discrepancy for set 2
for i = 1:size(v_Q2,1)
    % find closest point of Z2 to Q2
    [~,idx2] = min(abs(v_Z2-v_Q2(i)));
    d2(i) = f_Z2(idx2)-f_Q2(i);
end

% Maximum Vertical Deviation (MVD) and corresponding index
[MVD1,I1] = max(abs(d1));
[MVD2,I2] = max(abs(d2));

%% Plot first K-S test

% Define target axes
axes(handles.axes1)
% Plot ecdf of Z and Q (set1) 
hold on 
% Vertical line denotes maximum ecdf discrepancy
line([v_Q1(I1) v_Q1(I1)],[0 1],'Color',[1 1 0.5],'LineWidth',4)
% ecdf of Z1
stairs(v_Z1,f_Z1,'b','LineWidth',2)
% ecdf of Q1
stairs(v_Q1,f_Q1,'r','LineWidth',2)
axis([v_xmin v_xmax f_xmin f_xmax])
set(gca,'FontSize',9)
xlabel('File 1 & File 2 : Set 1','FontSize',8)
ylabel('Probability','FontSize',8)
grid on
box on 

% Set static text
set(handles.textKS1MVD,'String',num2str(MVD1,'%1.4f'))
set(handles.textKS1p,'String',num2str(p1,'%1.4f'))

%% Plot second K-S test

% Define target axes
axes(handles.axes2)
% Plot ecdf of Z and Q (set2)
hold on
% Vertical line denotes maximum ecdf discrepancy
line([v_Q2(I2) v_Q2(I2)],[0 1],'Color',[1 1 0.5],'LineWidth',4)
% ecdf of Z2
stairs(v_Z2,f_Z2,'b','LineWidth',2)
% ecdf of Q2
stairs(v_Q2,f_Q2,'r','LineWidth',2)
axis([v_ymin v_ymax f_ymin f_ymax])
set(gca,'FontSize',9)
xlabel('File 1 & File 2 : Set 2','FontSize',8)
ylabel('Probability','FontSize',8)
set(gca,'yaxislocation','right');
grid on
box on

% Set static text
set(handles.textKS2MVD,'String',num2str(MVD2,'%1.4f'))
set(handles.textKS2p,'String',num2str(p2,'%1.4f'))

%% CALCULATE SCATTERS (axis limits and correlations)

% compute max/min values for x and y axis
xmin = min([min(Z1) min(Q1)]);
xmax = max([max(Z1) max(Q1)]);
ymin = min([min(Z2) min(Q2)]);
ymax = max([max(Z2) max(Q2)]);

% adjust axis-limits for stationary variables
if xmin == xmax
    xmin = xmin - 1;
    xmax = xmax + 1;
end
if ymin == ymax
    ymin = ymin - 1;
    ymax = ymin + 1;
end

% Calculate corcondance measures
pearsonR1 = corr(Z1,Z2,'type','Pearson');
pearsonR2 = corr(Q1,Q2,'type','Pearson');
kendallTau1 = corr(Z1,Z2,'type','Kendall');
kendallTau2 = corr(Q1,Q2,'type','Kendall');
spearmanRho1 = corr(Z1,Z2,'type','Spearman');
spearmanRho2 = corr(Q1,Q2,'type','Spearman');

% Convex Hull computation
cv1 = convhull(Z1,Z2);
cv2 = convhull(Q1,Q2);

% Check which Q points lie within convex hull of Z
inConvexHull = inpolygon(Q1,Q2,Z1(cv1),Z2(cv1));

% Compute State-space coverage ratio (SSCR)
SSCR = sum(inConvexHull)/nObs2;

%%  PLOT SCATTERS

% Scatter plot historical data Z
axes(handles.axes3)
hold on
plot(Z1(cv1),Z2(cv1),'c')
plot(Q1(cv2),Q2(cv2),'m')
scatter(Z1,Z2,1,'b')
set(gca,'FontSize',9)
axis([xmin xmax ymin ymax])
grid on
box on
xlabel('File 1 : Set 1','FontSize',8)
ylabel('File 1 : Set 2','FontSize',8)

% Set static text
set(handles.textS1P,'String',num2str(pearsonR1,'%1.3f'))
set(handles.textS1K,'String',num2str(kendallTau1,'%1.3f'))
set(handles.textS1S,'String',num2str(spearmanRho1,'%1.3f'))

% Scatter Plot sampled data Q
axes(handles.axes4)
hold on
plot(Z1(cv1),Z2(cv1),'c')
plot(Q1(cv2),Q2(cv2),'m')
scatter(Q1,Q2,1,'r')
set(gca,'FontSize',9)
set(gca,'yaxislocation','right');
axis([xmin xmax ymin ymax])
grid on
box on  
xlabel('File 2 : Set 1','FontSize',8)
ylabel('File 2 : Set 2','FontSize',8)

% Set static text
set(handles.textS2P,'String',num2str(pearsonR2,'%1.3f'))
set(handles.textS2K,'String',num2str(kendallTau2,'%1.3f'))
set(handles.textS2S,'String',num2str(spearmanRho2,'%1.3f'))

%% CALCULATE HEATMAPS

% Heatmap will be displayed as an (HMsize x HMsize) square grid  
HMsize = 50; 

HM_Z = zeros(HMsize);
HM_Q = zeros(HMsize);

xstep = (xmax - xmin)/HMsize;
ystep = (ymax - ymin)/HMsize;

for i = 1:HMsize
    for j = 1:HMsize
        % define vertices of grid square
        gridVertices = [xmin + (i-1)*xstep  ymin + (j-1)*ystep;...
                        xmin + (i-1)*xstep  ymin +  j*ystep;...
                        xmin +  i*xstep     ymin +  j*ystep;...
                        xmin +  i*xstep     ymin + (j-1)*ystep];
        % compute heatmap of historical data Z
        HM_Z(j,i) = sum(inpolygon(Z1,Z2,gridVertices(:,1),gridVertices(:,2)));
        % Compute heatmap of sampled data Q
        HM_Q(j,i) = sum(inpolygon(Q1,Q2,gridVertices(:,1),gridVertices(:,2)));
    end
end

% Scale Heatmaps according to population
HM_Z = HM_Z/nObs1;
HM_Q = HM_Q/nObs2;

%%  PLOT HEATMAPS

% Heatmap of historical data Z
axes(handles.axes5)
imagesc(log10(HM_Z));
colormap('Gray');
set(gca,'YDir','normal') % flip Y axis back to normal
box on
% Remove XTicks and YTicks
set(gca,'XTick',[])
set(gca,'YTick',[])

% Heatmap of sampled data Q
axes(handles.axes6)
imagesc(log10(HM_Q));
colormap('Gray');
set(gca,'YDir','normal') % flip Y axis back to normal
box on
% Remove XTicks and YTicks
set(gca,'XTick',[])
set(gca,'YTick',[])


function fileName1TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to fileName1TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileName1TextBox as text
%        str2double(get(hObject,'String')) returns contents of fileName1TextBox as a double


% --- Executes during object creation, after setting all properties.
function fileName1TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileName1TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fileName2TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to fileName2TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileName2TextBox as text
%        str2double(get(hObject,'String')) returns contents of fileName2TextBox as a double



% --- Executes during object creation, after setting all properties.
function fileName2TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileName2TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browseButton1.
function browseButton1_Callback(hObject, eventdata, handles)
% hObject    handle to browseButton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName,pathName,~] = uigetfile('../*.mat');
fileName1 = fullfile(pathName,fileName);
set(handles.fileName1TextBox, 'String', fileName1)
Z = importdata(fileName1);
[nObs1,nVar1] = size(Z);
set(handles.nObs1TextBox, 'String', nObs1)
set(handles.nVar1TextBox, 'String', nVar1)

% --- Executes on button press in browseButton2.
function browseButton2_Callback(hObject, eventdata, handles)
% hObject    handle to browseButton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName,pathName,~] = uigetfile('../*.mat');
fileName2 = fullfile(pathName,fileName);
set(handles.fileName2TextBox, 'String', fileName2)
Q = importdata(fileName2);
[nObs2,nVar2] = size(Q);
set(handles.nObs2TextBox, 'String', nObs2)
set(handles.nVar2TextBox, 'String', nVar2)

function nObs1TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to nObs1TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nObs1TextBox as text
%        str2double(get(hObject,'String')) returns contents of nObs1TextBox as a double


% --- Executes during object creation, after setting all properties.
function nObs1TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nObs1TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nObs2TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to nObs2TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nObs2TextBox as text
%        str2double(get(hObject,'String')) returns contents of nObs2TextBox as a double


% --- Executes during object creation, after setting all properties.
function nObs2TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nObs2TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nVar1TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to nVar1TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nVar1TextBox as text
%        str2double(get(hObject,'String')) returns contents of nVar1TextBox as a double


% --- Executes during object creation, after setting all properties.
function nVar1TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nVar1TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nVar2TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to nVar2TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nVar2TextBox as text
%        str2double(get(hObject,'String')) returns contents of nVar2TextBox as a double


% --- Executes during object creation, after setting all properties.
function nVar2TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nVar2TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set1TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to set1TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set1TextBox as text
%        str2double(get(hObject,'String')) returns contents of set1TextBox as a double


% --- Executes during object creation, after setting all properties.
function set1TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set1TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set2TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to set2TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set2TextBox as text
%        str2double(get(hObject,'String')) returns contents of set2TextBox as a double


% --- Executes during object creation, after setting all properties.
function set2TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set2TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearButton.
function clearButton_Callback(hObject, eventdata, handles)
% hObject    handle to clearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clear all axes
cla(handles.axes1);
cla(handles.axes2);
cla(handles.axes3);
cla(handles.axes4);
cla(handles.axes5);
cla(handles.axes6);

% Clear static text
set(handles.textKS1MVD,'String','')
set(handles.textKS1p,'String','')

set(handles.textKS2MVD,'String','')
set(handles.textKS2p,'String','')


set(handles.textS1P,'String','')
set(handles.textS1K,'String','')
set(handles.textS1S,'String','')

set(handles.textS2P,'String','')
set(handles.textS2K,'String','')
set(handles.textS2S,'String','')

set(handles.set1StaticBox,'String','')
set(handles.set2StaticBox,'String','')

set(handles.nObs1TextBox,'String','')
set(handles.nVar1TextBox,'String','')
set(handles.nObs2TextBox,'String','')
set(handles.nVar2TextBox,'String','')


% --- Executes on button press in radiobuttonCustom.
function radiobuttonCustom_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonCustom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonCustom
set(handles.set1TextBox,'Enable','on')
set(handles.set2TextBox,'Enable','on')


% --- Executes on button press in radiobuttonDefault.
function radiobuttonDefault_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonDefault
set(handles.set1TextBox,'Enable','off')
set(handles.set2TextBox,'Enable','off')
