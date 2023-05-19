from distutils.cmd import Command
import tkinter
import customtkinter
from sys import exit
import os
from PIL import Image


#########################################################################################
#                              System Configurations
#########################################################################################
TITLE = "Quantum Satellite Simulator"
TITLE_FONT = ("Roboto", 45)
customtkinter.set_appearance_mode("dark")                                                 # Modes: system (default), light, dark
customtkinter.set_default_color_theme("blue")                                             # Themes: blue (default), dark-blue, green
#root = tkinter.Tk()
#SCREEN = [root.winfo_screenwidth(), root.winfo_screenheight()]

#########################################################################################
#                           Home Window Configurations
#########################################################################################
HOME_GEOMETRY = [640 , 162*4]                                                             # Window [width, height] in px
HOME_BTISPACE = 40                                                                        # Image to buttons spacing in px
HOME_BPADX    = 40                                                                        # Button X pading in px
HOME_BPADY    = 25                                                                        # Button Y pading in px
HOME_BSIZE    = [HOME_GEOMETRY[0]-2*HOME_BPADX , int(HOME_GEOMETRY[1]/4)-3*HOME_BPADY]    # Button [width, height] in px
HOME_BORDER   = 0                                                                         # Button border width in px
HOME_CORNER   = 8                                                                         # Button corner radius in px
HOME_FONT     = ("Roboto", 20)                                                            # Buttons Font
#########################################################################################
#                           Load Window Configurations
#########################################################################################
LOAD_GEOMETRY = [640 , 162*4]                                                             # Window [width, height] in px
#########################################################################################
#                         New Mission Window Configurations
#########################################################################################
NEW_TITLE = "New Project"
NEW_GEOMETRY = [690 , 640]                                                                # Window [width, height] in px 


class New_Mission(customtkinter.CTkToplevel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.title(NEW_TITLE)
        self.geometry('{}x{}'.format(NEW_GEOMETRY[0],NEW_GEOMETRY[0]))
        self.page = 0

        # set grid layout 4x4
        self.grid_rowconfigure(2, weight=1)
        self.grid_columnconfigure(4, weight=1)

        # Load Image
        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "GUI//Images")
        self.main_image = customtkinter.CTkImage(Image.open(image_path + "//Home_Qsat.jpg"),
                                               size=(HOME_GEOMETRY[0], int(HOME_GEOMETRY[0]/4)))
        # Insert Image
        self.home_image = customtkinter.CTkLabel(self, text="", image=self.main_image)
        self.home_image.grid(row=0, column=0, columnspan=4, padx=0, pady=0, sticky="nsew")

        # Create Frames
        self.first_frame()
        self.second_frame()
        self.third_frame()
        self.iframe.grid(row=1, column=0, columnspan=4, sticky="nsew")
        
        #Next and Back Button
        self.return_button = customtkinter.CTkButton(self,
                                                     text="Exit", command=self.back)
        self.return_button.grid(row=2, column=0, padx=20, pady=20, sticky="s")
        self.next_button = customtkinter.CTkButton(self,
                                                     text="Next", command=self.next)
        self.next_button.grid(row=2, column=3, padx=20, pady=20, sticky="s")

    def first_frame(self):
        self.iframe = customtkinter.CTkFrame(self, corner_radius=0, fg_color="transparent")
        self.iframe.grid_rowconfigure(5, weight=1)
        self.iframe.grid_columnconfigure(4, weight=1)

        # Name of the Mission
        self.nameLabel = customtkinter.CTkLabel(self.iframe, text="Mission Name")
        self.nameLabel.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
        self.nameEntry = customtkinter.CTkEntry(self.iframe, placeholder_text="Name")
        self.nameEntry.grid(row=0, column=1, columnspan=3, padx=20, pady=20, sticky="ew")
        # Number of states
        self.stateLabel = customtkinter.CTkLabel(self.iframe, text="Number of States")
        self.stateLabel.grid(row=1, column=0, padx=20, pady=20, sticky="nsew")
        self.stateEntry = customtkinter.CTkEntry(self.iframe, placeholder_text="4")
        self.stateEntry.grid(row=1, column=1, padx=20, pady=20, sticky="ew")
        # Number of decoys
        self.decoyLabel = customtkinter.CTkLabel(self.iframe, text="Number of Decoys")
        self.decoyLabel.grid(row=1, column=2, padx=20, pady=20, sticky="nsew")
        self.decoyEntry = customtkinter.CTkEntry(self.iframe, placeholder_text="2")
        self.decoyEntry.grid(row=1, column=3, padx=20, pady=20, sticky="ew")

        # Physical System Parameters
        self.phy_param_label = customtkinter.CTkLabel(self.iframe, text="Physical System Parameters",
                                                      font=customtkinter.CTkFont(size=20, weight="bold"))
        self.phy_param_label.grid(row=2, column=0, columnspan=4, padx=20, pady=20, sticky="nsew")
        #Orbit altitude h
        self.altLabel = customtkinter.CTkLabel(self.iframe, text="Orbit altitude")
        self.altLabel.grid(row=3, column=0, padx=20, pady=20, sticky="nsew")
        self.altEntry = customtkinter.CTkEntry(self.iframe, placeholder_text="m")
        self.altEntry.grid(row=3, column=1, padx=20, pady=20, sticky="ew")
        #Wavelenght wl
        self.altLabel = customtkinter.CTkLabel(self.iframe, text="Wavelenght")
        self.altLabel.grid(row=3, column=2, padx=20, pady=20, sticky="nsew")
        self.altEntry = customtkinter.CTkEntry(self.iframe, placeholder_text="nm")
        self.altEntry.grid(row=3, column=3, padx=20, pady=20, sticky="ew")
        #Transmitter aparture diameter dt
        self.altLabel = customtkinter.CTkLabel(self.iframe, text="Transmitter aparture")
        self.altLabel.grid(row=4, column=0, padx=20, pady=20, sticky="nsew")
        self.altEntry = customtkinter.CTkEntry(self.iframe, placeholder_text="cm")
        self.altEntry.grid(row=4, column=1, padx=20, pady=20, sticky="ew")
        #Receiver aparture diameter rt
        self.altLabel = customtkinter.CTkLabel(self.iframe, text="Receiver aparture")
        self.altLabel.grid(row=4, column=2, padx=20, pady=20, sticky="nsew")
        self.altEntry = customtkinter.CTkEntry(self.iframe, placeholder_text="cm")
        self.altEntry.grid(row=4, column=3, padx=20, pady=20, sticky="ew")
        #Gaussian Beam Waist w
        self.altLabel = customtkinter.CTkLabel(self.iframe, text="Beam Waist")
        self.altLabel.grid(row=5, column=0, padx=20, pady=20, sticky="nsew")
        self.altEntry = customtkinter.CTkEntry(self.iframe, placeholder_text="cm")
        self.altEntry.grid(row=5, column=1, padx=20, pady=20, sticky="ew")

    def second_frame(self):
        self.iiframe = customtkinter.CTkFrame(self, corner_radius=0, fg_color="transparent")
        self.iiframe.grid_rowconfigure(5, weight=1)
        self.iiframe.grid_columnconfigure(4, weight=1)

        self.mu1_Label = customtkinter.CTkLabel(self.iiframe, text="mu_1")
        self.mu1_Label.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
        self.mu1_Entry = customtkinter.CTkEntry(self.iiframe, placeholder_text="")
        self.mu1_Entry.grid(row=0, column=1, padx=20, pady=20, sticky="ew")
        self.mu1_Combo = customtkinter.CTkOptionMenu(self.iiframe, values=["Fixed", "Variable"],command=self.option_callback)
        self.mu1_Combo.grid(row=0, column=2, padx=20, pady=20, sticky="ew")
        self.mu1_Range = customtkinter.CTkEntry(self.iiframe, placeholder_text="<lower> , <upper>")

        self.mu2_Label = customtkinter.CTkLabel(self.iiframe, text="mu_2")
        self.mu2_Label.grid(row=1, column=0, padx=20, pady=20, sticky="nsew")
        self.mu2_Entry = customtkinter.CTkEntry(self.iiframe, placeholder_text="")
        self.mu2_Entry.grid(row=1, column=1, padx=20, pady=20, sticky="ew")
        self.mu2_Combo = customtkinter.CTkOptionMenu(self.iiframe, values=["Fixed", "Variable"],command=self.option_callback)
        self.mu2_Combo.grid(row=1, column=2, padx=20, pady=20, sticky="ew")
        self.mu2_Range = customtkinter.CTkEntry(self.iiframe, placeholder_text="<lower> , <upper>")

        self.mu3_Label = customtkinter.CTkLabel(self.iiframe, text="mu_3")
        self.mu3_Entry = customtkinter.CTkEntry(self.iiframe, placeholder_text="")
        self.mu3_Combo = customtkinter.CTkOptionMenu(self.iiframe, values=["Fixed", "Variable"],command=self.option_callback)
        self.mu3_Range = customtkinter.CTkEntry(self.iiframe, placeholder_text="<lower> , <upper>")

        self.pk1_Label = customtkinter.CTkLabel(self.iiframe, text="pk_1")
        self.pk1_Label.grid(row=3, column=0, padx=20, pady=20, sticky="nsew")
        self.pk1_Entry = customtkinter.CTkEntry(self.iiframe, placeholder_text="")
        self.pk1_Entry.grid(row=3, column=1, padx=20, pady=20, sticky="ew")
        self.pk1_Combo = customtkinter.CTkOptionMenu(self.iiframe, values=["Fixed", "Variable"],command=self.option_callback)
        self.pk1_Combo.grid(row=3, column=2, padx=20, pady=20, sticky="ew")
        self.pk1_Range = customtkinter.CTkEntry(self.iiframe, placeholder_text="<lower> , <upper>")

        self.pk2_Label = customtkinter.CTkLabel(self.iiframe, text="pk_2")
        self.pk2_Entry = customtkinter.CTkEntry(self.iiframe, placeholder_text="")
        self.pk2_Combo = customtkinter.CTkOptionMenu(self.iiframe, values=["Fixed", "Variable"],command=self.option_callback)
        self.pk2_Range = customtkinter.CTkEntry(self.iiframe, placeholder_text="<lower> , <upper>")

        self.pk2_Fill_Label = customtkinter.CTkLabel(self.iiframe, text="pk_2 = 1 - pk_1")

        self.pk3_Fill_Label = customtkinter.CTkLabel(self.iiframe, text="pk_3 = 1 - pk_2 - pk_1")

    def third_frame(self):
        self.iiiframe = customtkinter.CTkFrame(self, corner_radius=0, fg_color="transparent")
        self.iiiframe.grid_rowconfigure(5, weight=1)
        self.iiiframe.grid_columnconfigure(4, weight=1)

        self.pax_Label = customtkinter.CTkLabel(self.iiiframe, text="PaX")
        self.pax_Label.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
        self.pax_Entry = customtkinter.CTkEntry(self.iiiframe, placeholder_text="")
        self.pax_Entry.grid(row=0, column=1, padx=20, pady=20, sticky="ew")
        self.pbx_Label = customtkinter.CTkLabel(self.iiiframe, text="PbX")
        self.pbx_Label.grid(row=0, column=2, padx=20, pady=20, sticky="nsew")
        self.pbx_Entry = customtkinter.CTkEntry(self.iiiframe, placeholder_text="")
        self.pbx_Entry.grid(row=0, column=3, padx=20, pady=20, sticky="ew")

        self.ec_Label = customtkinter.CTkLabel(self.iiiframe, text="eps_c")
        self.ec_Label.grid(row=1, column=0, padx=20, pady=20, sticky="nsew")
        self.ec_Entry = customtkinter.CTkEntry(self.iiiframe, placeholder_text="")
        self.ec_Entry.grid(row=1, column=1, padx=20, pady=20, sticky="ew")
        self.es_Label = customtkinter.CTkLabel(self.iiiframe, text="eps_s")
        self.es_Label.grid(row=1, column=2, padx=20, pady=20, sticky="nsew")
        self.es_Entry = customtkinter.CTkEntry(self.iiiframe, placeholder_text="")
        self.es_Entry.grid(row=1, column=3, padx=20, pady=20, sticky="ew")

        self.pec_Label = customtkinter.CTkLabel(self.iiiframe, text="Pec")
        self.pec_Label.grid(row=2, column=0, padx=20, pady=20, sticky="nsew")
        self.pec_Entry = customtkinter.CTkEntry(self.iiiframe, placeholder_text="")
        self.pec_Entry.grid(row=2, column=1, padx=20, pady=20, sticky="ew")
        self.pap_Label = customtkinter.CTkLabel(self.iiiframe, text="Pap")
        self.pap_Label.grid(row=2, column=2, padx=20, pady=20, sticky="nsew")
        self.pap_Entry = customtkinter.CTkEntry(self.iiiframe, placeholder_text="")
        self.pap_Entry.grid(row=2, column=3, padx=20, pady=20, sticky="ew")

        self.Q_Label = customtkinter.CTkLabel(self.iiiframe, text="QBERI")
        self.Q_Label.grid(row=3, column=0, padx=20, pady=20, sticky="nsew")
        self.Q_Entry = customtkinter.CTkEntry(self.iiiframe, placeholder_text="")
        self.Q_Entry.grid(row=3, column=1, padx=20, pady=20, sticky="ew")
        self.srate_Label = customtkinter.CTkLabel(self.iiiframe, text="Source Rate")
        self.srate_Label.grid(row=3, column=2, padx=20, pady=20, sticky="nsew")
        self.srate_Entry = customtkinter.CTkEntry(self.iiiframe, placeholder_text="MHz")
        self.srate_Entry.grid(row=3, column=3, padx=20, pady=20, sticky="ew")

    def option_callback(self,choice):
        print(self.mu1_Combo.get())
        if self.mu1_Combo.get() == "Variable":
            self.mu1_Range.grid(row=0, column=3, padx=20, pady=20, sticky="ew")
        else:
            self.mu1_Range.grid_forget()
        if self.mu2_Combo.get() == "Variable":
            self.mu2_Range.grid(row=1, column=3, padx=20, pady=20, sticky="ew")
        else:
            self.mu2_Range.grid_forget()
        if self.mu3_Combo.get() == "Variable":
            self.mu3_Range.grid(row=2, column=3, padx=20, pady=20, sticky="ew")
        else:
            self.mu3_Range.grid_forget()
        if self.pk1_Combo.get() == "Variable":
            self.pk1_Range.grid(row=3, column=3, padx=20, pady=20, sticky="ew")
        else:
            self.pk1_Range.grid_forget()
        if self.pk2_Combo.get() == "Variable":
            self.pk2_Range.grid(row=4, column=3, padx=20, pady=20, sticky="ew")
        else:
            self.pk2_Range.grid_forget()

    def Verify_decoy(self):
        if int(self.decoyEntry.get()) > 1:
            self.pk2_Fill_Label.grid_forget()

            self.mu3_Label.grid(row=2, column=0, padx=20, pady=20, sticky="nsew")
            self.mu3_Entry.grid(row=2, column=1, padx=20, pady=20, sticky="ew")
            self.mu3_Combo.grid(row=2, column=2, padx=20, pady=20, sticky="ew")

            self.pk2_Label.grid(row=4, column=0, padx=20, pady=20, sticky="nsew")
            self.pk2_Entry.grid(row=4, column=1, padx=20, pady=20, sticky="ew")
            self.pk2_Combo.grid(row=4, column=2, padx=20, pady=20, sticky="ew")

            self.pk3_Fill_Label.grid(row=5, column=0, columnspan=2, padx=20, pady=20, sticky="nsew")
        else:
            self.mu3_Label.grid_forget()
            self.mu3_Entry.grid_forget()
            self.mu3_Combo.grid_forget()
            self.pk2_Label.grid_forget()
            self.pk2_Entry.grid_forget()
            self.pk2_Combo.grid_forget()
            self.pk3_Fill_Label_1.grid_forget()
            self.pk3_Fill_Label_2.grid_forget()

            self.pk2_Fill_Label.grid(row=4, column=0, columnspan=2, padx=20, pady=20, sticky="nsew")
            
    def select_frame_by_name(self, name):
        # show selected frame
        if name == "frame_1":
            self.iframe.grid(row=1, column=0, columnspan=4, sticky="nsew")
        else:
            self.iframe.grid_forget()
        if name == "frame_2":
            self.iiframe.grid(row=1, column=0, columnspan=4, sticky="nsew")
        else:
            self.iiframe.grid_forget()
        if name == "frame_3":
            self.iiiframe.grid(row=1, column=0, columnspan=4, sticky="nsew")
        else:
            self.iiiframe.grid_forget()

    def back(self):
        if self.page == 0:
            self.destroy()
        elif self.page == 1:
            self.page -= 1
            self.return_button.configure(text='Exit')
            self.select_frame_by_name('frame_1')
        elif self.page == 2:
            self.page -= 1
            self.next_button.configure(text="Next")
            self.select_frame_by_name('frame_2')
        return

    def next(self):
        if self.page == 0:
            self.page += 1
            self.return_button.configure(text='Back')
            self.Verify_decoy()
            self.select_frame_by_name('frame_2')
        elif self.page == 1:
            self.page += 1
            self.next_button.configure(text='Run Simulation')
            self.select_frame_by_name('frame_3')
        elif self.page == 2:
            #run sim
            return
        return

class Load_Mission(customtkinter.CTkToplevel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.title(TITLE)
        self.geometry('{}x{}'.format(LOAD_GEOMETRY[0],LOAD_GEOMETRY[0]))

        # set grid layout 1x2
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

        #Load Images
        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "GUI//Images")
        self.logo_image = customtkinter.CTkImage(Image.open(image_path + "//Logo.png"), size=(26, 26))

        # Load Folders
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Saves")
        folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]

        # create navigation frame
        self.navigation_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.navigation_frame.grid(row=0, column=0, sticky="nsew")
        self.navigation_frame.grid_rowconfigure(len(folders)+1, weight=1)

        self.navigation_frame_label = customtkinter.CTkLabel(self.navigation_frame, text="  Saved Missions", image=self.logo_image,
                                                             compound="left", font=customtkinter.CTkFont(size=15, weight="bold"))
        self.navigation_frame_label.grid(row=0, column=0, padx=20, pady=20)

        self.file_button = []
        for i in range(len(folders)):
            self.file_button += [customtkinter.CTkButton(self.navigation_frame, corner_radius=0, height=40, border_spacing=10, text=str(folders[i]),
                                                       fg_color="transparent", text_color=("gray10", "gray90"), hover_color=("gray70", "gray30"),
                                                       anchor="w", command=self.select(folders[i]))]
            self.file_button[i].grid(row=i+1, column=0, sticky="ew")

        self.return_button = customtkinter.CTkButton(self.navigation_frame,
                                                     text="Back", command=self.back)
        self.return_button.grid(row=6, column=0, padx=20, pady=20, sticky="s")

    def select(self,folder):
        return

    def back(self):
        self.destroy()
        return

class Home(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.title(TITLE)
        self.geometry('{}x{}'.format(HOME_GEOMETRY[0],HOME_GEOMETRY[0]))

        # Load Image
        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "GUI//Images")
        self.main_image = customtkinter.CTkImage(Image.open(image_path + "//Home_Qsat.jpg"),
                                               size=(HOME_GEOMETRY[0], int(HOME_GEOMETRY[0]/4)))

        # Insert Image
        self.home_image = customtkinter.CTkLabel(self, text="", image=self.main_image)
        self.home_image.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")

        # Use Buttons
        self.button_New = customtkinter.CTkButton(master=self, width=HOME_BSIZE[0], height=HOME_BSIZE[1],
                                                  border_width=HOME_BORDER, corner_radius=HOME_CORNER,
                                                  text="New Mission", font=HOME_FONT, command=self.b_new)
        self.button_New.grid(row=1, column=0, padx=HOME_BPADX, pady=(HOME_BPADY+50,HOME_BPADY), sticky="ew")
        self.button_Load = customtkinter.CTkButton(master=self, width=HOME_BSIZE[0], height=HOME_BSIZE[1],
                                                  border_width=HOME_BORDER, corner_radius=HOME_CORNER,
                                                  text="Load Mission", font=HOME_FONT, command=self.b_load)
        self.button_Load.grid(row=2, column=0, padx=HOME_BPADX, pady=HOME_BPADY, sticky="ew")
        self.button_Exit = customtkinter.CTkButton(master=self, width=HOME_BSIZE[0], height=HOME_BSIZE[1],
                                                  border_width=HOME_BORDER, corner_radius=HOME_CORNER,
                                                  text="Exit", font=HOME_FONT, command=self.b_exit)
        self.button_Exit.grid(row=3, column=0, padx=HOME_BPADX, pady=HOME_BPADY, sticky="ew")

        # Windows
        self.new_mission = None
        self.load_mission = None


    def b_new(self):
        if self.new_mission is None or not self.new_mission.winfo_exists():
            self.new_mission = New_Mission(self)  # create window if its None or destroyed
        else:
            self.new_mission.focus()  # if window exists focus it
        return

    def b_load(self):
        if self.load_mission is None or not self.load_mission.winfo_exists():
            self.load_mission = Load_Mission(self)  # create window if its None or destroyed
        else:
            self.load_mission.focus()  # if window exists focus it
        return

    def b_exit(self):
        exit(1)

if __name__ == "__main__":
    app = Home()
    app.mainloop()

