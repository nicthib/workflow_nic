function varargout = VideoSitichGUI(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VideoSitichGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @VideoSitichGUI_OutputFcn, ...
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

function VideoSitichGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.MainDir = '/local_mount/space/revault/revault2/cmdata_CCD_analysis/cm62_2';
handles.nWebcam = 0;
handles.nData = 0;
guidata(hObject, handles);

function varargout = VideoSitichGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



function LoadWebcamVid_Callback(hObject, eventdata, handles)
data.{['webcam' mat2str(handles.nWebcam+1)]} = load(uigetfile(handles.MainDir),'webcam');
handles.nWebcam = handles.nWebcam + 1;

guidata(hObject, handles);




function LoadDataVid_Callback(hObject, eventdata, handles)

function PreviewBtn_Callback(hObject, eventdata, handles)

function MakeVIdeo_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Crop1_CreateFcn(hObject, eventdata, handles)

function Crop2_CreateFcn(hObject, eventdata, handles)

function Crop3_CreateFcn(hObject, eventdata, handles)

function Crop4_CreateFcn(hObject, eventdata, handles)

function canvasW_CreateFcn(hObject, eventdata, handles)

function canvasH_CreateFcn(hObject, eventdata, handles)


