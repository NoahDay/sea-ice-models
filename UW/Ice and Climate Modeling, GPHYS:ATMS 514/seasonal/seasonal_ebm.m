function fig = seasonal_ebm()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

load seasonal_ebm
global guihnd;

h0 = figure('Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'FileName','ebm.fig.m', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'Position',[19 527 825 397], ...
	'Tag','Fig1', ...
	'ToolBar','none');
guihnd=h0;

h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[404.3992987204725 348.7719365157481 59.5552411417323 20.17193651574803], ...
	'String','D =', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat1, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[472.5996555118111 342.2176205708661 99.84 30.72], ...
	'String','0.44', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat2, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[671.04 341.76 99.84 30.72], ...
	'String','1', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[602.88 348.48 59.52 20.16], ...
	'String','Q/Qo =', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat3, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[472.5996555118111 287.44 99.89911417322836 29.77762057086614], ...
	'String','203.3', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[404.3992987204725 293.6 59.5552411417323 20.17193651574803], ...
	'String','A =', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat4, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[472.5996555118111 232.44 100 30], ...
	'String','2.09', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[404.3992987204725 238.44 59.52 20.16], ...
	'String','B = ', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback',mat5, ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[233.4181225393701 72.04263041338584 141.2035556102362 36.50159940944882], ...
	'String','use defaults', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback',mat6, ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[46.10728346456693 70.12149360236221 161.3754921259843 87.41172490157481], ...
	'String','Run EBM', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback',mat7, ...
	'FontSize',12, ...
	'ListboxTop',0, ...
	'Position',[233.4181225393701 289.92 137.28 32.64], ...
	'String','No Albedo Feedback', ...
	'Style','checkbox', ...
	'Tag','Checkbox1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback',mat8, ...
	'FontSize',12, ...
	'ListboxTop',0, ...
	'Position',[233.4181225393701 249.6 137.28 34.56], ...
	'String','Simulate Hadley Cell', ...
	'Style','checkbox', ...
	'Tag','Checkbox1', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback',mat9, ...
	'FontSize',12, ...
	'ListboxTop',0, ...
	'Position',[233.4181225393701 213.12 137.28 31.68], ...
	'String','Use Cold Start', ...
	'Style','checkbox', ...
	'Tag','Checkbox1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[602.88 183.36 59.52 20.16], ...
	'String','per =', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat10, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[671.04 177.6 99.84 29.76], ...
	'String','102.07', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[602.88 238.08 59.52 20.16], ...
	'String','obl =', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat11, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[671.04 232.32 99.84 29.76], ...
	'String','23.44', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat12, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[671.04 287.04 99.84 29.76], ...
	'String','0.01672', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[602.88 293.76 59.52 20.16], ...
	'String','ecc =', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'FontSize',12, ...
	'ListboxTop',0, ...
	'Position',[37.44 180.48 179.52 143.04], ...
	'String','Period and Land Distribution', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat13, ...
	'FontSize',12, ...
	'Position',[41.28000000000002 194.88 171.84 102.72], ...
	'String',mat14, ...
	'Style','listbox', ...
	'Tag','Listbox1', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback',mat15, ...
	'FontSize',12, ...
	'ListboxTop',0, ...
	'Position',[233.4181225393701 176.64 137.28 31.68], ...
	'String','Explicit Sea Ice', ...
	'Style','checkbox', ...
	'Tag','Checkbox1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[37.44 341.76 98.88 22.08], ...
	'String','Case Name =', ...
	'Style','text', ...
	'Tag','StaticText3');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat16, ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[146.88 338.88 226.56 26.88], ...
	'String','Control', ...
	'Style','edit', ...
	'Tag','EditText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat17, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',mat18, ...
	'String','30', ...
	'Style','edit', ...
	'Tag','EditText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',mat19, ...
	'String','Run Length =', ...
	'Style','text', ...
	'Tag','StaticText3');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',mat20, ...
	'String','Gridpoints =', ...
	'Style','text', ...
	'Tag','StaticText3');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat21, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[721.3868725393702 74.92433562992127 46.10728346456693 26.89591535433071], ...
	'String','50', ...
	'Style','edit', ...
	'Tag','EditText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[404.3992987204725 72.95999999999999 59.52 20.16], ...
	'String','Cw =', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat22, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[472.5996555118111 66.23999999999999 99.84 29.76], ...
	'String','9.8', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[404.3992987204725 128.12 59.52 20.16], ...
	'String','Cl =', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat23, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[472.5996555118111 121 99.84 30.72], ...
	'String','0.45', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat24, ...
	'FontSize',18, ...
	'ListboxTop',0, ...
	'Position',[472.5996555118111 176.72 99.84 30.72], ...
	'String','3', ...
	'Style','edit', ...
	'Tag','EditText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'FontSize',14, ...
	'ListboxTop',0, ...
	'Position',[404.3992987204725 183.28 59.52 20.16], ...
	'String','nu =', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.701960784313725 0.701960784313725 0.701960784313725], ...
	'Callback',mat25, ...
	'FontSize',12, ...
	'ListboxTop',0, ...
	'Position',[233.4181225393701 121.9921875 140.2429872047244 34.5804625984252], ...
	'String','Save annual cycle', ...
	'Tag','Pushbutton2');

if nargout > 0, fig = h0; end
