import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import HeatLib

class HeatEquationApp:
    PARAMETER_LABELS = [
        "Initial Temperature",
        "Left Temperature parameter a:",
        "Right Temperature parameter a:",
        "Left Heat Flux f:",
        "Right Heat Flux f:",
        "Left Heat Coefficient k:",
        "Right Heat Coefficient k:",
        "Spatial step dx:",
        "Time step dt:",
        "L:",
        "T:"
    ]

    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Heat Equation App")
        self.root.geometry("500x300")
        self.method_var = tk.StringVar(value="FTCS")
        self.entries = []
        self.create_main_window()

    def create_main_window(self):
        main_frame = tk.Frame(self.root)
        main_frame.pack(anchor="w", padx=20, pady=20, fill="x")

        tk.Label(main_frame,text="Select Solving Method:",font=("Helvetica", 14)).pack(anchor="w")

        for method_name in ["FTCS", "Crank-Nicolson"]:
            tk.Radiobutton(main_frame,text=method_name,variable=self.method_var,value=method_name,font=("Helvetica", 12)).pack(anchor="w")

        tk.Button(main_frame,text="Submit",command=self.submit_method).pack(anchor="w", pady=10)

    def submit_method(self):
        method = self.method_var.get()
        self.root.withdraw()
        self.open_parameter_window(method)

    def open_parameter_window(self, method):
        bg_color = "#f0f0f0"

        self.parameter_window = tk.Toplevel(self.root)
        self.parameter_window.title(f"{method} Parameters")
        self.parameter_window.geometry("1000x650")
        self.parameter_window.config(bg=bg_color)
        self.parameter_window.columnconfigure(1, weight=1)

        self.entries = []

        fig, ax = plt.subplots(figsize=(4, 1))
        fig.patch.set_facecolor(bg_color)
        ax.set_facecolor(bg_color)
        ax.axis("off")

        equation = r"$k\frac{\partial u(x_i,t)}{\partial x}=f-au_i$"
        ax.text(0.5, 0.5, equation, fontsize=28, ha="center", va="center")

        canvas = FigureCanvasTkAgg(fig, master=self.parameter_window)
        canvas.get_tk_widget().grid(row=0, column=1, pady=10)

        tk.Label(self.parameter_window,text="Boundary Conditions",font=("Helvetica", 18),bg=bg_color).grid(row=0, column=0, padx=20, sticky="w")

        for row, label_text in enumerate(self.PARAMETER_LABELS, start=1):

            tk.Label(self.parameter_window,text=label_text,font=("Helvetica", 16),bg=bg_color).grid(row=row, column=0, sticky="w", padx=20, pady=2)

            entry = tk.Entry(self.parameter_window,font=("Helvetica", 16),width=40)
            entry.grid(row=row, column=1, sticky="w", pady=2)

            self.entries.append(entry)

        tk.Button(self.parameter_window,text="Run",command=lambda: self.run_solver(method)).grid(row=13, column=1, sticky="w", pady=15)

        self.parameter_window.protocol("WM_DELETE_WINDOW", self.root.destroy)

    def run_solver(self, method):
        try:
            input_values = (
                self.entries[0].get(),
                float(self.entries[1].get()),
                float(self.entries[2].get()),
                self.entries[3].get(),
                self.entries[4].get(),
                float(self.entries[5].get()),
                float(self.entries[6].get()),
                float(self.entries[7].get()),
                float(self.entries[8].get()),
                float(self.entries[9].get()),
                float(self.entries[10].get())
            )

            if method == "FTCS":
                HeatLib.FTCS(*input_values)
            elif method == "Crank-Nicolson":
                HeatLib.CrankNicolson(*input_values)

            self.parameter_window.destroy()

        except ValueError:
            messagebox.showwarning("Input Error","Please enter valid numeric values.")
        except Exception as error:messagebox.showerror("Error",f"Something went wrong:\n{error}")

    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    app = HeatEquationApp()
    app.run()
