import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.font_manager as fm
from matplotlib.figure import Figure

class CatppuccinTheme:
    # Catppuccin Mocha palette
    ROSEWATER = "#f5e0dc"
    FLAMINGO = "#f2cdcd"
    PINK = "#f5c2e7"
    MAUVE = "#cba6f7"
    RED = "#f38ba8"
    MAROON = "#eba0ac"
    PEACH = "#fab387"
    YELLOW = "#f9e2af"
    GREEN = "#a6e3a1"
    TEAL = "#94e2d5"
    SKY = "#89dceb"
    SAPPHIRE = "#74c7ec"
    BLUE = "#89b4fa"
    LAVENDER = "#b4befe"
    TEXT = "#cdd6f4"
    SUBTEXT1 = "#bac2de"
    SUBTEXT0 = "#a6adc8"
    OVERLAY2 = "#9399b2"
    OVERLAY1 = "#7f849c"
    OVERLAY0 = "#6c7086"
    SURFACE2 = "#585b70"
    SURFACE1 = "#45475a"
    SURFACE0 = "#313244"
    BASE = "#1e1e2e"
    MANTLE = "#181825"
    CRUST = "#11111b"

class ModernODEAnalyzer:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("ODE Analysis Suite")
        self.root.configure(bg=CatppuccinTheme.BASE)
        
        # Configure styles
        self.setup_styles()
        
        # Create main frame
        self.frame = ttk.Frame(self.root, padding="20", style='Main.TFrame')
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Title
        title_label = ttk.Label(
            self.frame,
            text="ODE Analysis Suite",
            style='Title.TLabel'
        )
        title_label.grid(row=0, column=0, columnspan=2, pady=(0, 20))
        
        # ODE Type Selection
        self.ode_type = tk.StringVar(value="1d_autonomous")
        ttk.Label(
            self.frame,
            text="Select ODE Type:",
            style='Header.TLabel'
        ).grid(row=1, column=0, pady=5)
        
        self.type_menu = ttk.Combobox(
            self.frame,
            textvariable=self.ode_type,
            style='Modern.TCombobox'
        )
        self.type_menu['values'] = ('1d_autonomous', '2nd_order', '2d_system')
        self.type_menu.grid(row=1, column=1, pady=5, padx=5, sticky='ew')
        self.type_menu.bind('<<ComboboxSelected>>', self.update_parameters)
        
        # Equation Display
        self.equation_frame = ttk.LabelFrame(
            self.frame,
            text="Current Equation",
            padding="10",
            style='Equation.TLabelframe'
        )
        self.equation_frame.grid(row=2, column=0, columnspan=2, pady=10, sticky='ew')
        self.equation_label = ttk.Label(
            self.equation_frame,
            text="",
            style='Equation.TLabel'
        )
        self.equation_label.pack(pady=5)
        
        # Parameters frame
        self.param_frame = ttk.LabelFrame(
            self.frame,
            text="Parameters",
            padding="10",
            style='Modern.TLabelframe'
        )
        self.param_frame.grid(row=3, column=0, columnspan=2, pady=10, sticky='ew')
        
        # Create figure with Catppuccin theme
        plt.style.use('dark_background')
        self.fig = Figure(figsize=(8, 5), facecolor=CatppuccinTheme.BASE)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor(CatppuccinTheme.SURFACE0)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.get_tk_widget().grid(row=4, column=0, columnspan=2, pady=10)
        
        # Initial parameter setup
        self.param_vars = {}
        self.setup_1d_autonomous()
        
        # Solve button with modern styling
        solve_button = ttk.Button(
            self.frame,
            text="Solve",
            command=self.solve,
            style='Accent.TButton'
        )
        solve_button.grid(row=5, column=0, columnspan=2, pady=10)

    def setup_styles(self):
        style = ttk.Style()
        
        # Configure modern styles
        style.configure('Main.TFrame',
                       background=CatppuccinTheme.BASE)
        
        style.configure('Title.TLabel',
                       font=('Helvetica', 24, 'bold'),
                       foreground=CatppuccinTheme.LAVENDER,
                       background=CatppuccinTheme.BASE)
        
        style.configure('Header.TLabel',
                       font=('Helvetica', 12),
                       foreground=CatppuccinTheme.TEXT,
                       background=CatppuccinTheme.BASE)
        
        style.configure('Equation.TLabel',
                       font=('Courier', 12),
                       foreground=CatppuccinTheme.GREEN,
                       background=CatppuccinTheme.SURFACE0)
        
        style.configure('Equation.TLabelframe',
                       background=CatppuccinTheme.SURFACE0,
                       foreground=CatppuccinTheme.TEXT)
        
        style.configure('Modern.TLabelframe',
                       background=CatppuccinTheme.SURFACE0,
                       foreground=CatppuccinTheme.TEXT)
        
        style.configure('Modern.TLabelframe.Label',
                       font=('Helvetica', 10, 'bold'),
                       foreground=CatppuccinTheme.TEXT,
                       background=CatppuccinTheme.SURFACE0)
        
        style.configure('Accent.TButton',
                       font=('Helvetica', 12),
                       padding=10,
                       background=CatppuccinTheme.MAUVE,
                       foreground=CatppuccinTheme.BASE)
        
        # Make entry widgets round and modern
        style.configure('Modern.TEntry',
                       fieldbackground=CatppuccinTheme.SURFACE1,
                       foreground=CatppuccinTheme.TEXT,
                       padding=5)

    def update_equation_display(self):
        equations = {
            '1d_autonomous': "dx/dt = rx(1 - x/K)",
            '2nd_order': "mx'' + cx' + kx = 0",
            '2d_system': "dx/dt = a₁₁x + a₁₂y\ndy/dt = a₂₁x + a₂₂y"
        }
        self.equation_label.config(text=equations[self.ode_type.get()])

    def create_modern_entry(self, parent, variable, row, column, label_text):
        ttk.Label(
            parent,
            text=label_text,
            style='Header.TLabel'
        ).grid(row=row, column=column, padx=5, pady=5)
        
        entry = ttk.Entry(
            parent,
            textvariable=variable,
            style='Modern.TEntry'
        )
        entry.grid(row=row, column=column+1, padx=5, pady=5, sticky='ew')
        return entry

    def setup_1d_autonomous(self):
        for widget in self.param_frame.winfo_children():
            widget.destroy()
        
        self.param_vars = {
            'r': tk.DoubleVar(value=2.0),
            'K': tk.DoubleVar(value=100.0),
            'x0': tk.DoubleVar(value=10.0),
            'dt': tk.DoubleVar(value=0.1),
            'T': tk.DoubleVar(value=10.0)
        }
        
        row = 0
        for name, var in self.param_vars.items():
            self.create_modern_entry(self.param_frame, var, row, 0, f"{name}:")
            row += 1
            
        self.update_equation_display()

    def setup_2nd_order(self):
        for widget in self.param_frame.winfo_children():
            widget.destroy()
            
        self.param_vars = {
            'mass': tk.DoubleVar(value=1.0),
            'spring_const': tk.DoubleVar(value=10.0),
            'damping': tk.DoubleVar(value=0.5),
            'x0': tk.DoubleVar(value=1.0),
            'v0': tk.DoubleVar(value=0.0),
            'T': tk.DoubleVar(value=10.0)
        }
        
        row = 0
        for name, var in self.param_vars.items():
            self.create_modern_entry(self.param_frame, var, row, 0, f"{name}:")
            row += 1
            
        self.update_equation_display()

    def setup_2d_system(self):
        for widget in self.param_frame.winfo_children():
            widget.destroy()
            
        self.param_vars = {
            'a11': tk.DoubleVar(value=2.0),
            'a12': tk.DoubleVar(value=-1.0),
            'a21': tk.DoubleVar(value=1.0),
            'a22': tk.DoubleVar(value=-2.0),
            'x0': tk.DoubleVar(value=1.0),
            'y0': tk.DoubleVar(value=0.0),
            'T': tk.DoubleVar(value=10.0)
        }
        
        row = 0
        for name, var in self.param_vars.items():
            self.create_modern_entry(self.param_frame, var, row, 0, f"{name}:")
            row += 1
            
        self.update_equation_display()

    def update_parameters(self, event=None):
        if self.ode_type.get() == '1d_autonomous':
            self.setup_1d_autonomous()
        elif self.ode_type.get() == '2nd_order':
            self.setup_2nd_order()
        else:
            self.setup_2d_system()

    def plot_solution(self, t, x, y=None):
        self.ax.clear()
        self.ax.set_facecolor(CatppuccinTheme.SURFACE0)
        
        if y is None:
            self.ax.plot(t, x, color=CatppuccinTheme.MAUVE, linewidth=2)
        else:
            self.ax.plot(t, x, color=CatppuccinTheme.MAUVE, linewidth=2, label='x(t)')
            self.ax.plot(t, y, color=CatppuccinTheme.GREEN, linewidth=2, label='y(t)')
            self.ax.legend()
        
        self.ax.set_xlabel('Time', color=CatppuccinTheme.TEXT)
        self.ax.set_ylabel('Value', color=CatppuccinTheme.TEXT)
        self.ax.grid(True, color=CatppuccinTheme.SURFACE1)
        self.ax.tick_params(colors=CatppuccinTheme.TEXT)
        
        # Update spines
        for spine in self.ax.spines.values():
            spine.set_color(CatppuccinTheme.SURFACE1)
        
        self.canvas.draw()

    def solve_1d_autonomous(self):
        r = self.param_vars['r'].get()
        K = self.param_vars['K'].get()
        x0 = self.param_vars['x0'].get()
        dt = self.param_vars['dt'].get()
        T = self.param_vars['T'].get()
        
        t = np.arange(0, T, dt)
        x = np.zeros_like(t)
        x[0] = x0
        
        for i in range(1, len(t)):
            dxdt = r * x[i-1] * (1 - x[i-1]/K)
            x[i] = x[i-1] + dt * dxdt
        
        self.plot_solution(t, x)

    def solve_2nd_order(self):
        m = self.param_vars['mass'].get()
        k = self.param_vars['spring_const'].get()
        c = self.param_vars['damping'].get()
        x0 = self.param_vars['x0'].get()
        v0 = self.param_vars['v0'].get()
        T = self.param_vars['T'].get()
        
        omega = np.sqrt(k/m)
        zeta = c/(2*np.sqrt(m*k))
        
        t = np.linspace(0, T, 1000)
        if zeta < 1:  # Underdamped
            omega_d = omega * np.sqrt(1 - zeta**2)
            x = np.exp(-zeta*omega*t) * (x0*np.cos(omega_d*t) + 
                (v0 + zeta*omega*x0)*np.sin(omega_d*t)/omega_d)
        elif zeta == 1:  # Critically damped
            x = x0*(1 + omega*t)*np.exp(-omega*t)
        else:  # Overdamped
            r1 = -omega*(zeta + np.sqrt(zeta**2 - 1))
            r2 = -omega*(zeta - np.sqrt(zeta**2 - 1))
            x = (x0*(r1*np.exp(r2*t) - r2*np.exp(r1*t)))/(r1 - r2)
        
        self.plot_solution(t, x)

    def solve_2d_system(self):
        A = np.array([[self.param_vars['a11'].get(), self.param_vars['a12'].get()],
                     [self.param_vars['a21'].get(), self.param_vars['a22'].get()]])
        x0 = np.array([self.param_vars['x0'].get(), self.param_vars['y0'].get()])
        T = self.param_vars['T'].get()
        
        eigenvals, eigenvecs = np.linalg.eig(A)
        
        t = np.linspace(0, T, 1000)
        X = np.zeros((len(t), 2))
        
        if np.isreal(eigenvals).all():  # Real eigenvalues
            c = np.linalg.solve(eigenvecs, x0)
            for i, ti in enumerate(t):
                X[i] = c[0]*eigenvecs[:,0]*np.exp(eigenvals[0]*ti) + \
                    c[1]*eigenvecs[:,1]*np.exp(eigenvals[1]*ti)

        else:  # Complex eigenvalues
            a = eigenvals[0].real
            b = eigenvals[0].imag
            for i, ti in enumerate(t):
                X[i] = np.exp(a*ti) * (x0*np.cos(b*ti) + (A@x0 - a*x0)/b*np.sin(b*ti))
                
        self.plot_solution(t, X[:,0], X[:,1])

    def solve(self):
        if self.ode_type.get() == '1d_autonomous':
            self.solve_1d_autonomous()
        elif self.ode_type.get() == '2nd_order':
            self.solve_2nd_order()
        else:
            self.solve_2d_system()

    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    app = ModernODEAnalyzer()
    app.run()
