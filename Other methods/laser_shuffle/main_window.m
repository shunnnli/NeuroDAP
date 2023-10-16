function varargout = main_window(varargin)
% MAIN_WINDOW MATLAB code for main_window.fig
%      MAIN_WINDOW, by itself, creates a new MAIN_WINDOW or raises the existing
%      singleton*.
%
%      H = MAIN_WINDOW returns the handle to a new MAIN_WINDOW or the handle to
%      the existing singleton*.
%
%      MAIN_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_WINDOW.M with the given input arguments.
%
%      MAIN_WINDOW('Property','Value',...) creates a new MAIN_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_window

% Last Modified by GUIDE v2.5 02-Oct-2012 12:56:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_window_OpeningFcn, ...
                   'gui_OutputFcn',  @main_window_OutputFcn, ...
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

controller = LaserShuffleController;

% --- Executes just before main_window is made visible.
function main_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_window (see VARARGIN)

% Choose default command line output for main_window
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_data_dir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_data_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_data_dir as text
%        str2double(get(hObject,'String')) returns contents of edit_data_dir as a double


% --- Executes during object creation, after setting all properties.
function edit_data_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_data_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_browse_data_dir.
function button_browse_data_dir_Callback(hObject, eventdata, handles)
% hObject    handle to button_browse_data_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    controller.dataDir = uigetdir();
    


% --- Executes on selection change in listbox_valid_files.
function listbox_valid_files_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_valid_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_valid_files contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_valid_files


% --- Executes during object creation, after setting all properties.
function listbox_valid_files_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_valid_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_laser_name.
function listbox_laser_name_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_laser_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_laser_name contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_laser_name


% --- Executes during object creation, after setting all properties.
function listbox_laser_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_laser_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
