import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import HeatLib

class HeatEquationApp:

    PARAMETER_LABELS = [
        "Initial Temperature:",
        "Left parameter a:",
        "Right parameter a:",
        "Initial Left Boundary Value (n=0):",
        "Initial Right Boundary Value (n=0):",
        "New Left Boundary f (n>=1):",
        "New Right Boundary f (n>=1):",
        "Left coefficient k:",
        "Right coefficient k:",
        "Spatial step dx:",
        "Time step dt:",
        "Domain length L:",
        "Simulation time T:"
    ]

# ----------------------------------------------------------------------------------------------------------------------------------------

    def __init__(self):

        self.root = tk.Tk()
        self.root.title("Heat Equation App")
        self.root.geometry("500x300")

        self.method_var = tk.StringVar(value="FTCS")
        self.entries = []

        self.create_main_window()

# ----------------------------------------------------------------------------------------------------------------------------------------

    def create_main_window(self):

        main_frame = tk.Frame(self.root)
        main_frame.pack(anchor="w", padx=20, pady=20, fill="x")

        tk.Label(main_frame,text="Select Solving Method:",font=("Helvetica", 14)).pack(anchor="w")

        for method_name in ["FTCS", "Crank-Nicolson"]:

            tk.Radiobutton(main_frame,text=method_name,variable=self.method_var,value=method_name,font=("Helvetica", 12)).pack(anchor="w")

        tk.Button(main_frame,text="Submit",command=self.submit_method).pack(anchor="w", pady=10)

# ------------------------------------------------------------------------------------------------------------------------------------------

    def submit_method(self):

        method = self.method_var.get()
        self.root.withdraw()
        self.open_parameter_window(method)

# ------------------------------------------------------------------------------------------------------------------------------------------

    def open_parameter_window(self, method):

        bg_color = "#f0f0f0"

        self.parameter_window = tk.Toplevel(self.root)
        self.parameter_window.title(f"{method} Parameters")
        self.parameter_window.geometry("1200x680")
        self.parameter_window.minsize(1050, 600)
        self.parameter_window.configure(bg=bg_color)

        # --- Main two-column container -------------------------------

        main_frame = tk.Frame(self.parameter_window,bg=bg_color)

        main_frame.pack(fill="both",expand=True,padx=20,pady=20)

        main_frame.columnconfigure(0, weight=3)
        main_frame.columnconfigure(1, weight=2)

        # --- LEFT: Parameters ----------------------------------------

        left_frame = tk.Frame(main_frame,bg=bg_color)

        left_frame.grid(row=0,column=0,sticky="nsew",padx=(0, 25))

        tk.Label(
            left_frame,
            text="Simulation Parameters",
            font=("Helvetica", 18, "bold"),
            bg=bg_color
        ).grid(
            row=0,
            column=0,
            columnspan=2,
            sticky="w",
            pady=(0, 15)
        )

        self.entries = []

        for row, label_text in enumerate(self.PARAMETER_LABELS,start=1):

            tk.Label(
                left_frame,
                text=label_text,
                font=("Helvetica", 11),
                bg=bg_color
            ).grid(
                row=row,
                column=0,
                sticky="w",
                pady=5,
                padx=(0, 12)
            )

            entry = tk.Entry(left_frame,font=("Helvetica", 11),width=25)

            entry.grid(row=row,column=1,sticky="ew",pady=5)

            self.entries.append(entry)

        left_frame.columnconfigure(1, weight=1)

        # --- Run button ----------------------------------------------

        tk.Button(
            left_frame,
            text="Run Simulation",
            font=("Helvetica", 12, "bold"),
            command=lambda: self.run_solver(method),
            width=20
        ).grid(
            row=14,
            column=1,
            sticky="w",
            pady=(20, 0)
        )

        # --- RIGHT: Boundary information -----------------------------

        right_frame = tk.Frame(main_frame,bg="white",bd=1,relief="solid")

        right_frame.grid(row=0,column=1,sticky="nsew")

        tk.Label(
            right_frame,
            text="Boundary Conditions",
            font=("Helvetica", 17, "bold"),
            bg="white"
        ).pack(
            anchor="w",
            padx=20,
            pady=(20, 5)
        )

        # --- Mathematical equation ----------------------------------

        fig, ax = plt.subplots(figsize=(4.5, 1.1))

        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")
        ax.axis("off")

        equation = (r"$a_ u + k" r"\frac{\partial u}{\partial n}=f$")

        ax.text(0.5,0.5, equation,fontsize=21,ha="center",va="center")

        canvas = FigureCanvasTkAgg(fig, master=right_frame)

        canvas.get_tk_widget().pack(fill="x", padx=10)

        plt.close(fig)

        # --- Boundary condition guide -------------------------------

        guide = (
            "DIRICHLET  (k = 0)\n"
            "Initial: a × initial boundary temperature\n"
            "New: a × prescribed boundary temperature\n\n"

            "NEUMANN  (a = 0)\n"
            "Initial: k × initial outward-normal gradient\n"
            "New: k × prescribed outward-normal gradient\n\n"

            "ROBIN  (a != 0, k != 0)\n"
            "Initial: a × initial temperature\n"
            "             + k × initial outward-normal gradient\n"
            "New: prescribed Robin forcing f(t)\n\n"

            "OUTWARD-NORMAL CONVENTION\n"
            "Left boundary:  -du/dx\n"
            "Right boundary: +du/dx\n\n"

            "IMPORTANT\n"
            "Calculate the initial boundary values from\n"
            "your initial temperature profile before\n"
            "running the simulation.\n\n"
            "HeatLib preserves the initial temperature\n"
            "distribution at n = 0."
        )

        tk.Label(
            right_frame,
            text=guide,
            font=("Helvetica", 11),
            bg="white",
            justify="left",
            anchor="nw",
            wraplength=410
        ).pack(
            anchor="w",
            fill="x",
            padx=20,
            pady=10
        )

        self.parameter_window.protocol("WM_DELETE_WINDOW",self.root.destroy)

# ---------------------------------------------------------------------------------------------------------------------------------

    def run_solver(self, method):

        try:

            values = [entry.get().strip() for entry in self.entries]

            TempFormula = values[0]

            a0 = float(values[1])
            aN = float(values[2])

            f0_initial = float(values[3])
            fN_initial = float(values[4])

            f0 = values[5]
            fN = values[6]

            k0 = float(values[7])
            kN = float(values[8])

            dx = float(values[9])
            dt = float(values[10])
            L = float(values[11])
            T = float(values[12])

            if method == "FTCS":

                HeatLib.FTCS(TempFormula, a0, aN, f0_initial, fN_initial, f0, fN, k0, kN, dx, dt, L, T)

            elif method == "Crank-Nicolson":

                # Uses the existing CN interface until updated.
                HeatLib.CrankNicolson(TempFormula,a0, aN,f0, fN, k0, kN, dx, dt, L, T)

            self.parameter_window.destroy()

        except ValueError as error:

            messagebox.showwarning("Input Error",str(error))

        except Exception as error:

            messagebox.showerror("Error", f"Something went wrong:\n{error}")

# ---------------------------------------------------------------------------------------------------------------------------------------

    def run(self):

        self.root.mainloop()


if __name__ == "__main__":

    app = HeatEquationApp()
    app.run()

