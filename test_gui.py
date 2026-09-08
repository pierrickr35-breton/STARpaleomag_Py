import tkinter as tk
from tkinter import ttk, messagebox

class STARpaleomagApp:
    def __init__(self, root):
        self.root = root
        self.root.title("STARpaleomag_Py - Paléomagnétisme")
        self.root.geometry("1100x700")

        # 1. Barre de menus supérieure
        self._setup_menu()

        # 2. Panneau séparateur réajustable (Gauche: Graphique, Droite: Texte)
        self.paned_window = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        self.paned_window.pack(fill=tk.BOTH, expand=True)

        # --- PANNEAU GAUCHE : Zone Graphique (Rendu Canvas / SVG) ---
        self.graph_frame = ttk.Frame(self.paned_window, width=550)
        self.canvas = tk.Canvas(self.graph_frame, bg="white")
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self.paned_window.add(self.graph_frame, weight=1)

        # --- PANNEAU DROIT : Console / Liste de données ---
        self.text_frame = ttk.Frame(self.paned_window, width=550)
        self.text_area = tk.Text(
            self.text_frame, 
            bg="#1b1c1e", 
            fg="#ffffff", 
            insertbackground="white", 
            font=("Courier", 10)
        )
        self.text_area.pack(fill=tk.BOTH, expand=True)
        self.paned_window.add(self.text_frame, weight=1)

        # Initialisation des tracés et des données
        self.draw_demo_plots()
        self.load_demo_data()

    def _setup_menu(self):
        menubar = tk.Menu(self.root)

        # Menu Fichiers
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Ouvrir...", command=self.dummy_action)
        file_menu.add_separator()
        file_menu.add_command(label="Quitter", command=self.root.quit)
        menubar.add_cascade(label="Fichiers", menu=file_menu)

        # Menu Selection données
        data_menu = tk.Menu(menubar, tearoff=0)
        data_menu.add_command(label="Sélectionner échantillon", command=self.dummy_action)
        menubar.add_cascade(label="Selection donnees", menu=data_menu)

        # Menu Graphiques (Reproduisant la capture d'écran)
        graph_menu = tk.Menu(menubar, tearoff=0)
        graph_menu.add_command(label="Zijderveld", command=self.draw_demo_plots)
        graph_menu.add_command(label="Paramètres Zijderveld", command=self.dummy_action)
        graph_menu.add_separator()
        graph_menu.add_command(label="XYGraph", command=self.dummy_action)
        graph_menu.add_separator()
        graph_menu.add_command(label="Stereo Mesures", command=self.dummy_action)
        graph_menu.add_command(label="Stereo Résultats", command=self.dummy_action)
        graph_menu.add_command(label="Paramètres Stereo", command=self.dummy_action)
        graph_menu.add_separator()
        graph_menu.add_command(label="Susceptibilité", command=self.dummy_action)
        graph_menu.add_command(label="Données+interprétation", command=self.dummy_action)
        graph_menu.add_separator()
        graph_menu.add_command(label="Clear Screen", command=self.clear_screen)
        menubar.add_cascade(label="Graphiques", menu=graph_menu)

        # Menu Calcul & Graphics_SVG
        menubar.add_cascade(label="Calcul", menu=tk.Menu(menubar, tearoff=0))
        menubar.add_cascade(label="Graphics_SVG", menu=tk.Menu(menubar, tearoff=0))

        self.root.config(menu=menubar)

    def draw_demo_plots(self):
        """Exemple d'affichage graphique (Stéréogramme + Zijderveld)"""
        self.clear_screen()
        
        # 1. Stéréogramme (Haut)
        self.canvas.create_oval(150, 20, 350, 220, outline="black", width=2)
        self.canvas.create_line(250, 20, 250, 220, fill="gray", dash=(2, 2))
        self.canvas.create_line(150, 120, 350, 120, fill="gray", dash=(2, 2))
        self.canvas.create_text(250, 10, text="N", font=("Arial", 10, "bold"))
        self.canvas.create_text(250, 230, text="S", font=("Arial", 10, "bold"))
        self.canvas.create_text(140, 120, text="W", font=("Arial", 10, "bold"))
        self.canvas.create_text(360, 120, text="E", font=("Arial", 10, "bold"))
        
        # Points sur le stéréogramme
        for (x, y) in [(290, 130), (292, 132), (295, 129)]:
            self.canvas.create_oval(x-3, y-3, x+3, y+3, outline="black", fill="white")

        # 2. Diagramme de Zijderveld (Bas)
        self.canvas.create_line(250, 280, 250, 480, width=2) # Axe vertical (Down)
        self.canvas.create_line(100, 430, 400, 430, width=2) # Axe horizontal (E)
        self.canvas.create_text(250, 495, text="Down", font=("Arial", 10, "bold"))
        self.canvas.create_text(415, 430, text="E", font=("Arial", 10, "bold"))
        
        # Projection horizontale (Points verts)
        pts_green = [(120, 410), (150, 370), (180, 330), (210, 290), (230, 260)]
        for i in range(len(pts_green)-1):
            self.canvas.create_line(pts_green[i], pts_green[i+1], fill="green", width=2)
            self.canvas.create_oval(pts_green[i][0]-3, pts_green[i][1]-3, pts_green[i][0]+3, pts_green[i][1]+3, outline="green")

        # Projection verticale (Points rouges)
        pts_red = [(120, 430), (160, 440), (200, 450), (240, 460), (280, 465)]
        for i in range(len(pts_red)-1):
            self.canvas.create_line(pts_red[i], pts_red[i+1], fill="red", width=2)
            self.canvas.create_oval(pts_red[i][0]-3, pts_red[i][1]-3, pts_red[i][0]+3, pts_red[i][1]+3, outline="red", fill="red")

    def load_demo_data(self):
        """Simulation de la console de texte avec les données Pmag"""
        header = "List dataPmag\n" + " Numero".ljust(15) + "Dec".rjust(8) + "Inc".rjust(8) + "q".rjust(6) + "Mag".rjust(8) + "\n"
        header += "-" * 50 + "\n"
        self.text_area.insert(tk.END, header)

        data = [
            ("1: AP2-8a", "108.8", "-49.5", "0", "nd"),
            ("2: AP2-8a", "108.4", "-49.9", "0", "nd"),
            ("3: AP2-8a", "109.2", "-50.5", "0", "nd"),
            ("4: AP2-8a", "109.1", "-53.2", "0", "nd"),
            ("5: AP2-8a", "113.4", "-55.0", "0", "nd"),
            ("6: AP2-8a", "110.3", "-54.9", "0", "nd"),
        ]

        for num, dec, inc, q, mag in data:
            line = f"{num.ljust(15)}{dec.rjust(8)}{inc.rjust(8)}{q.rjust(6)}{mag.rjust(8)}\n"
            self.text_area.insert(tk.END, line)

    def clear_screen(self):
        self.canvas.delete("all")

    def dummy_action(self):
        pass

if __name__ == "__main__":
    root = tk.Tk()
    app = STARpaleomagApp(root)
    root.mainloop()