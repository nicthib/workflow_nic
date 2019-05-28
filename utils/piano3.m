% FRATRES: (C#, A, E), (Bb, E, D), (A, E, C#), (G, C, Bb), (F, C, A), (E, C, G), (D, C, F) and (C#, A, E)
function varargout = piano3(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @piano3_OpeningFcn, ...
                   'gui_OutputFcn',  @piano3_OutputFcn, ...
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

function piano3_OpeningFcn(hObject, eventdata, handles, varargin)
addpath('/Users/Nic/Desktop/MIDI code/Code');
handles.keys = []; handles.allnotes = []; handles.nchord = 1;
handles.thisfullchord = [];
handles.practiceflag = 0;
handles.nnotes = 18;
handles.output = hObject;
handles.octOffset = 3;
set(handles.OctText,'string',mat2str(handles.octOffset));
guidata(hObject, handles);

function varargout = piano3_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function C_Callback(hObject, eventdata, handles)
handles = Presskey(handles,0);
guidata(hObject,handles)

function CSharp_Callback(hObject, eventdata, handles)
handles = Presskey(handles,1);
guidata(hObject,handles)

function D_Callback(hObject, eventdata, handles)
handles = Presskey(handles,2);
guidata(hObject,handles)

function DSharp_Callback(hObject, eventdata, handles)
handles = Presskey(handles,3);
guidata(hObject,handles)

function E_Callback(hObject, eventdata, handles)
handles = Presskey(handles,4);
guidata(hObject,handles)

function F_Callback(hObject, eventdata, handles)
handles = Presskey(handles,5);
guidata(hObject,handles)

function FSharp_Callback(hObject, eventdata, handles)
handles = Presskey(handles,6);
guidata(hObject,handles)

function G_Callback(hObject, eventdata, handles)
handles = Presskey(handles,7);
guidata(hObject,handles)

function GSharp_Callback(hObject, eventdata, handles)
handles = Presskey(handles,8);
guidata(hObject,handles)

function A_Callback(hObject, eventdata, handles)
handles = Presskey(handles,9);
guidata(hObject,handles)

function ASharp_Callback(hObject, eventdata, handles)
handles = Presskey(handles,10);
guidata(hObject,handles)

function B_Callback(hObject, eventdata, handles)
handles = Presskey(handles,11);
guidata(hObject,handles)

function C2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,12);
guidata(hObject,handles)

function CSharp2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,13);
guidata(hObject,handles)

function D2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,14);
guidata(hObject,handles)

function DSharp2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,15);
guidata(hObject,handles)

function E2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,16);
guidata(hObject,handles)

function F2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,17);
guidata(hObject,handles)

function FSharp2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,18);
guidata(hObject,handles)

function G2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,19);
guidata(hObject,handles)

function GSharp2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,20);
guidata(hObject,handles)

function A2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,21);
guidata(hObject,handles)

function ASharp2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,22);
guidata(hObject,handles)

function B2_Callback(hObject, eventdata, handles)
handles = Presskey(handles,23);
guidata(hObject,handles)


function edit4_Callback(hObject, eventdata, handles)
handles.nnotes = str2double(get(hObject,'String'));
handles.keys = [];
guidata(hObject,handles)

function nchords_Callback(hObject, eventdata, handles)
handles.nChords = [];
handles.nChords = str2double(get(hObject,'String'));
handles.allnotes = zeros(handles.nChords,50);
handles.keys = [];

guidata(hObject,handles)

function OctUp_Callback(hObject, eventdata, handles)
handles.octOffset = handles.octOffset + 1;
set(handles.OctText,'string',mat2str(handles.octOffset));
guidata(hObject,handles)

function OctDown_Callback(hObject, eventdata, handles)
handles.octOffset = handles.octOffset - 1;
set(handles.OctText,'string',mat2str(handles.octOffset));
guidata(hObject,handles)

function NextChord_Callback(hObject, eventdata, handles)
handles.thisfullchord = [];
for i = 1:ceil(handles.nnotes/length(handles.keys))
   handles.thisfullchord = [handles.thisfullchord handles.keys+12*(i-1)];
end
handles.allnotes(handles.nchord,1:length(handles.thisfullchord)) = handles.thisfullchord;
handles.nchord = handles.nchord + 1;

if handles.nchord>handles.nChords
    close_Callback(hObject, eventdata, handles)
    return
end
handles.keys = [];
guidata(hObject,handles)

function ResetChord_Callback(hObject, eventdata, handles)
handles.thisfullchord = [];
guidata(hObject,handles)

function close_Callback(hObject, eventdata, handles)
assignin('base','keys',handles.allnotes);
close all;

function handles = Presskey(handles,Note)
if handles.practiceflag == 0
	handles.keys = [handles.keys handles.octOffset*12 + Note];
end
PlayNote(handles.octOffset*12 + Note)

function PlayNote(n)
f = midi2freq(n);
Fs = 44100;
y = synth(f,.5,.3,Fs,'sine');
sound(y,Fs);

function radiobutton6_Callback(hObject, eventdata, handles)
handles.practiceflag = get(hObject,'Value');
guidata(hObject,handles)

function nchords_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ChordStr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
