import matplotlib.pyplot as plt
import numpy as np
import sys

# Default Plot Settings:
DEFAULT_SET = {}
DEFAULT_SET['style'] = 'seaborn-colorblind'
DEFAULT_SET['size'] = (10,6)
DEFAULT_SET['linewidth'] = 3
DEFAULT_SET['linestyle'] = ''
DEFAULT_SET['font'] = 13
DEFAULT_SET['label location'] = 'upper right'
DEFAULT_SET['grid'] = False
DEFAULT_SET['twincolor'] = ['blue', 'red']

class plot:

    def __init__(self, settings = None):
        
        if settings is not None and isinstance(settings, plt_set):
            self.default_settings = settings
        else:
            self.default_settings = plt_set()

        self.data = []
        self.twindata = []

        self.ttl = text( "Title" , self.default_settings.font )
        self.x_leg = text( "X Legend" , self.default_settings.font )
        self.y_leg = text( "Y Legend" , self.default_settings.font )
        self.yaux_leg = text( "Y Legend" , self.default_settings.font )

    def title(self, title_text, font=None):
        self.ttl = text( title_text , self.default_settings.font )
        if font is not None:
            self.ttl.Font(font)

    def x_legend(self, legend_text, font=None):
        self.x_leg = text( legend_text , self.default_settings.font )
        if font is not None:
            self.x_leg.Font(font)

    def y_legend(self, legend_text, font=None):
        self.y_leg = text( legend_text , self.default_settings.font )
        if font is not None:
            self.y_leg.Font(font)

    def ytwin_legend(self, legend_text, font=None):
        self.yaux_leg = text( legend_text , self.default_settings.font )
        if font is not None:
            self.yaux_leg.Font(font)

    def add_data(self, x_values, y_values, label=None, invert=None, style=None, linewidth=None):
        self.data += [ data_set(x_values, y_values, None, [False,False],
                                self.default_settings.linewidth, self.default_settings.linestyle) ]
        
        if style is not None:
            self.data[-1].Line_style(style)
        if linewidth is not None:
            self.data[-1].Line_width(linewidth)
        if label is not None:
            self.data[-1].Label(label)
        if invert is not None:
            self.data[-1].Invert(invert)

    def add_twindata(self, x_values, y_values, label=None, invert=None, style=None, linewidth=None):
        self.twindata += [ data_set(x_values, y_values, None, [False,False],
                                    self.default_settings.linewidth, self.default_settings.linestyle) ]
        if style is not None:
            self.twindata[-1].Line_style(style)
        if linewidth is not None:
            self.twindata[-1].Line_width(linewidth)
        if label is not None:
            self.twindata[-1].Label(label)
        if invert is not None:
            self.twindata[-1].Invert(invert)

    def pop_data(self):
        self.data.pop()
    
    def pop_all(self):
        self.data = []

    def pop_twindata(self):
        self.twindata.pop()

    def plot(self):
        
        plt.style.use(self.default_settings.style)
        plt.rcParams["figure.figsize"] = self.default_settings.size

        for i in range(len(self.data)):
            plt.plot(self.data[i].x, self.data[i].y, self.data[i].linestyle,
                     linewidth=self.data[i].linewidth, label=self.data[i].label)

        plt.xlabel(self.x_leg.text, fontsize=self.x_leg.font)
        plt.ylabel(self.y_leg.text, fontsize=self.y_leg.font)
        plt.title(self.ttl.text, fontsize=self.ttl.font)

        if ifnone([x.label for x in self.data]):
            plt.legend(loc=self.default_settings.loc)

        if self.default_settings.grid:
            plt.grid()

        if self.data[0].invert[0]:
            plt.invert_xaxis()

        if self.data[0].invert[1]:
            plt.invert_yaxis()
        
        plt.show()

    def logplot(self):

        plt.style.use(self.default_settings.style)
        plt.rcParams["figure.figsize"] = self.default_settings.size

        for i in range(len(self.data)):
            plt.plot(self.data[i].x, self.data[i].y, self.data[i].linestyle,
                     linewidth=self.data[i].linewidth, label=self.data[i].label)

        plt.xlabel(self.x_leg.text, fontsize=self.x_leg.font)
        plt.ylabel(self.y_leg.text, fontsize=self.y_leg.font)
        plt.yscale("log")
        plt.title(self.ttl.text, fontsize=self.ttl.font)

        if ifnone([x.label for x in self.data]):
            plt.legend(loc=self.default_settings.loc)

        if self.default_settings.grid:
            plt.grid()

        if self.data[0].invert[0]:
            plt.invert_xaxis()

        if self.data[0].invert[1]:
            plt.invert_yaxis()
        
        plt.show()

    def twinplot(self, scale=None):
        # Create figure and axis objects
        fig, ax1 = plt.subplots()

        # Plot Style
        plt.style.use(self.default_settings.style)
        plt.rcParams["figure.figsize"] = self.default_settings.size

        # Plot the first data set on the left y-axis
        for i in range(len(self.data)):
            ax1.plot(self.data[i].x, self.data[i].y, self.data[i].linestyle,
                     linewidth=self.data[i].linewidth, label=self.data[i].label)
        
        ax1.set_xlabel(self.x_leg.text, fontsize=self.x_leg.font)
        ax1.set_ylabel(self.y_leg.text, fontsize=self.y_leg.font, color=self.default_settings.twincolor[0])
        ax1.tick_params(axis='y', labelcolor=self.default_settings.twincolor[0])

        # Create a second y-axis with a different scale
        ax2 = ax1.twinx()
        for i in range(len(self.twindata)):
            ax2.plot(self.twindata[i].x, self.twindata[i].y, self.twindata[i].linestyle,
                     linewidth=self.twindata[i].linewidth, label=self.twindata[i].label)
        
        ax2.set_ylabel(self.yaux_leg.text, fontsize=self.yaux_leg.font, color=self.default_settings.twincolor[1])
        ax2.tick_params(axis='y', labelcolor=self.default_settings.twincolor[1])
        
        if scale is not None:
            ax1.set_yscale('log')
        
        if self.default_settings.grid:
            ax1.grid()

        if self.data[0].invert[0]:
            plt.invert_xaxis()

        if self.data[0].invert[1]:
            ax1.invert_yaxis()

        if self.twindata[0].invert[1]:
            ax2.invert_yaxis()

        # Add title and legend
        plt.title(self.ttl.text, fontsize=self.ttl.font)
        if ifnone([x.label for x in self.data]):
            ax1.legend(loc="upper left")
        if ifnone([x.label for x in self.twindata]):
            ax2.legend(loc="upper right")

        # Show the plot
        plt.show()

class plt_set:

    def __init__(self):
        self.style = DEFAULT_SET['style']
        self.size = DEFAULT_SET['size']
        self.font = DEFAULT_SET['font']
        self.linewidth = DEFAULT_SET['linewidth']
        self.linestyle = DEFAULT_SET['linestyle']
        self.loc = DEFAULT_SET['label location']
        self.grid = DEFAULT_SET['grid']

        self.twincolor = DEFAULT_SET['twincolor']

    def Style(self,s):
        self.style = s

    def Size(self,s):
        if s is tuple:
            self.size = s

    def Font(self,f):
        self.font = f

    def LineWidth(self,lw):
        self.linewidth = lw

    def LineStyle(self,ls):
        self.linestyle = ls

    def Grid(self,g):
        self.grid = g

    def Clear(self):
        self.__init__()  

class data_set:

    def __init__(self,x_values,y_values,label,invert,linewidth,linestyle):
        
        if dim(x_values)[0] == dim(y_values)[0]:

            self.x = np.array(x_values)
            self.y = np.array(y_values)
            self.linewidth = linewidth
            self.linestyle = linestyle
            self.label = label
            self.invert = invert
        
        else:
            sys.stderr.write( "Error on data set: X ({}) and Y ({}) have diffrent dimentions".format(dim(x_values)[0],dim(y_values)[0]) )

    def Line_width(self,lw):
        self.linewidth = lw

    def Line_style(self,ls):
        self.linestyle = ls

    def Label(self,l):
        self.label = l

    def Invert(self,I):
        self.invert = I

class text:

    def __init__(self, text , font):
        self.text = text
        self.font = font

    def Font(self,f):
        self.font = f


def dim(arr):
    if isinstance(arr, list):
        if len(arr) > 0:
            return [len(arr)] + dim(arr[0])
        else:
            return []
    elif isinstance(arr, np.ndarray):
        return list(arr.shape)
    else:
        return []

def ifnone(arr):

    i=0
    flag = False
    while i < len(arr) and flag == False:
        if arr[i] is not None:
            flag = True
        i += 1

    return flag


