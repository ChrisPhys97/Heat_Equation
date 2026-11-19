import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import HeatLib


class HeatEquationApp:

    PARAMETER_LABELS = [
        "Initial Temperature:",
        "Left Temperature parameter a:",
        "Right Temperature parameter a:",
        "Initial Left Boundary Value (n=0):",
        "Initial Right Boundary Value (n=0):",
        "New Left Boundary f (n>=1):",
        "New Right Boundary f (n>=1):",
        "Left Heat Coefficient k:",
        "Right Heat Coefficient k:",
        "Spatial step dx:",
        "Time step dt:",
        "L:",
        "T:"
    ]

 # ---------------------------------------------------------------------
    
    def __init__(self):

        self.root = tk.Tk()
        self.root.title("Heat Equation App")
        self.root.geometry("500x300")

        self.method_var = tk.StringVar(value="FTCS")
        self.entries = []

        self.create_main_window()

# ------------------------------------------------------------------------

    def create_main_window(self):

        main_frame = tk.Frame(self.root)
        main_frame.pack(anchor="w",padx=20,pady=20,fill="x")

        tk.Label(main_frame,text="Select Solving Method:",font=("Helvetica", 14)).pack(anchor="w")

        for method_name in ["FTCS", "Crank-Nicolson"]:

            tk.Radiobutton(
                main_frame,
                text=method_name,
                variable=self.method_var,
                value=method_name,
                font=("Helvetica", 12)
            ).pack(anchor="w")

        tk.Button(main_frame,text="Submit",command=self.submit_method).pack(anchor="w", pady=10)

# -------------------------------------------------------------------------------------------------------------

    def submit_method(self):

        method = self.method_var.get()

        self.root.withdraw()
        self.open_parameter_window(method)

# --------------------------------------------------------------------------------------------------------------

    def open_parameter_window(self, method):

        bg_color = "#f0f0f0"

        self.parameter_window = tk.Toplevel(self.root)
        self.parameter_window.title(f"{method} Parameters")
        self.parameter_window.geometry("1050x850")
        self.parameter_window.config(bg=bg_color)

        self.parameter_window.columnconfigure(1, weight=1)

        self.entries = []

        # --- Display general boundary equation ------------------------

        fig, ax = plt.subplots(figsize=(5, 1))

        fig.patch.set_facecolor(bg_color)
        ax.set_facecolor(bg_color)
        ax.axis("off")

        equation = (r"$a_b u + k_b" r"\frac{\partial u}{\partial n}=f_b$")

        ax.text(0.5, 0.5,equation,fontsize=24,ha="center",va="center")

        canvas = FigureCanvasTkAgg(fig,master=self.parameter_window)

        canvas.get_tk_widget().grid(row=0,column=1,pady=10)

        # --- Heading -------------------------------------------------

        tk.Label(
            self.parameter_window,
            text="Boundary Conditions",
            font=("Helvetica", 18),
            bg=bg_color
        ).grid(
            row=0,
            column=0,
            padx=20,
            sticky="w"
        )

        # --- Parameter entries ---------------------------------------

        for row, label_text in enumerate(
            self.PARAMETER_LABELS,
            start=1
        ):

            tk.Label(
                self.parameter_window,
                text=label_text,
                font=("Helvetica", 14),
                bg=bg_color
            ).grid(
                row=row,
                column=0,
                sticky="w",
                padx=20,
                pady=3
            )

            entry = tk.Entry(
                self.parameter_window,
                font=("Helvetica", 14),
                width=35
            )

            entry.grid(
                row=row,
                column=1,
                sticky="w",
                pady=3
            )

            self.entries.append(entry)

        # --- Boundary condition explanation --------------------------

        explanation = (
            "BOUNDARY CONDITION GUIDE\n\n"
            "Dirichlet (k = 0):\n"
            "  Initial value: a × initial boundary temperature.\n"
            "  New value: a × prescribed boundary temperature.\n\n"
            "Neumann (a = 0):\n"
            "  Initial value: k × initial outward-normal gradient.\n"
            "  New value: k × prescribed outward-normal gradient.\n\n"
            "Robin (a != 0, k != 0):\n"
            "  Initial value: a × initial temperature\n"
            "                 + k × initial outward-normal gradient.\n"
            "  New value: prescribed Robin forcing f(t).\n\n"
            "IMPORTANT:\n"
            "Calculate the initial boundary values from your initial\n"
            "temperature profile before running the simulation.\n"
            "Left outward gradient = -du/dx\n"
            "Right outward gradient = +du/dx\n"
            "The initial temperature profile is not modified."
        )

        tk.Label(
            self.parameter_window,
            text=explanation,
            font=("Helvetica", 10),
            bg=bg_color,
            justify="left",
            anchor="w"
        ).grid(
            row=14,
            column=0,
            columnspan=2,
            sticky="w",
            padx=20,
            pady=10
        )

        # --- Run button ----------------------------------------------

        tk.Button(
            self.parameter_window,
            text="Run",
            font=("Helvetica", 14),
            command=lambda: self.run_solver(method)
        ).grid(
            row=15,
            column=1,
            sticky="w",
            pady=15
        )

        self.parameter_window.protocol(
            "WM_DELETE_WINDOW",
            self.root.destroy
        )

    def run_solver(self, method):

        try:

            values = [entry.get().strip() for entry in self.entries]

            # Initial temperature formula
            TempFormula = values[0]

            # Boundary coefficients
            a0 = float(values[1])
            aN = float(values[2])

            # Initial boundary forcing values
            f0_initial = float(values[3])
            fN_initial = float(values[4])

            # New boundary forcing expressions
            f0 = values[5]
            fN = values[6]

            k0 = float(values[7])
            kN = float(values[8])

            # Grid and simulation parameters
            dx = float(values[9])
            dt = float(values[10])
            L = float(values[11])
            T = float(values[12])

            if method == "FTCS":

                HeatLib.FTCS(TempFormula, a0, aN, f0_initial, fN_initial, f0, fN, k0, kN, dx, dt, L, T)

            elif method == "Crank-Nicolson":

                # Existing CN signature:
                # TempFormula, a0, aN, f0, fN,
                # k0, kN, dx, dt, L, T

                HeatLib.CrankNicolson(TempFormula, a0, aN, f0, fN, k0, kN, dx, dt, L, T)

            self.parameter_window.destroy()

        except ValueError as error:

            messagebox.showwarning("Input Error", str(error))

        except Exception as error:

            messagebox.showerror("Error",f"Something went wrong:\n{error}")

    def run(self):
        self.root.mainloop()

if __name__ == "__main__":

    app = HeatEquationApp()
    app.run()

