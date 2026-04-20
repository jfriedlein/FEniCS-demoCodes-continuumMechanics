import os
import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import subprocess
import threading
import re
import sys
import platform

# Image support
try:
    from PIL import Image, ImageTk
    PIL_AVAILABLE = True
except:
    PIL_AVAILABLE = False

# 🔥 Match run.sh root
BASE_DIR = os.getcwd()

# Detect GUI
try:
    root = tk.Tk()
    root.withdraw()
    GUI_MODE = True
except tk.TclError:
    GUI_MODE = False
    print("01-elasto-static/01-tractions")
    sys.exit(0)


# ---------------- Image Loader ---------------- #

def load_image(path, size=(80, 80)):
    if not PIL_AVAILABLE:
        return None
    try:
        img = Image.open(path)
        img.thumbnail(size)
        return ImageTk.PhotoImage(img)
    except:
        return None


# ---------------- MAIN CLASS ---------------- #

class ExampleSelector:
    def __init__(self, root):
        self.root = root
        self.root.title("FEniCS Example Selector")
        self.root.geometry("1000x720")
        self.root.configure(bg="#f0f0f0")
        self.root.deiconify()

        # 🎨 Style
        style = ttk.Style()
        style.theme_use("clam")

        style.configure("Card.TFrame",
                        background="white",
                        relief="raised",
                        borderwidth=1)

        style.configure("TButton",
                        font=("Arial", 10),
                        padding=6)

        self.main_frame = ttk.Frame(root, padding=20)
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        self.show_topics()

    # ---------------- Helpers ---------------- #

    def clear(self):
        for w in self.main_frame.winfo_children():
            w.destroy()

    def safe_print(self, text):
        self.root.after(0, lambda: self.output.insert(tk.END, text))

    # ---------------- Topics ---------------- #

    def show_topics(self):
        self.clear()

        ttk.Label(self.main_frame, text="Select Topic",
                  font=("Arial", 20, "bold")).pack(pady=10)

        frame = ttk.Frame(self.main_frame)
        frame.pack(fill=tk.BOTH, expand=True)

        topics = [
            d for d in os.listdir(BASE_DIR)
            if os.path.isdir(os.path.join(BASE_DIR, d)) and re.match(r'^\d{2}-', d)
        ]
        topics.sort()

        cols = 3
        for i, topic in enumerate(topics):
            r, c = divmod(i, cols)

            card = ttk.Frame(frame, style="Card.TFrame", padding=10)
            card.grid(row=r, column=c, padx=15, pady=15, sticky="nsew")

            img = load_image(os.path.join(BASE_DIR, topic, "icon.png"))
            if img:
                lbl = ttk.Label(card, image=img)
                lbl.image = img
                lbl.pack()

            ttk.Button(card, text=topic,
                       command=lambda t=topic: self.open_topic(t)).pack(fill=tk.X)

        for c in range(cols):
            frame.grid_columnconfigure(c, weight=1)

    # ---------------- Smart Topic ---------------- #

    def open_topic(self, topic):
        topic_path = os.path.join(BASE_DIR, topic)

        if os.path.exists(os.path.join(topic_path, "main.py")):
            self.show_run(topic, None)
        else:
            self.show_examples(topic)

    # ---------------- Examples ---------------- #

    def show_examples(self, topic):
        self.clear()

        ttk.Button(self.main_frame, text="← Back",
                   command=self.show_topics).pack(anchor="w")

        ttk.Label(self.main_frame, text=topic,
                  font=("Arial", 18, "bold")).pack(pady=10)

        frame = ttk.Frame(self.main_frame)
        frame.pack(fill=tk.BOTH, expand=True)

        topic_path = os.path.join(BASE_DIR, topic)

        examples = [
            d for d in os.listdir(topic_path)
            if os.path.isdir(os.path.join(topic_path, d))
        ]
        examples.sort()

        cols = 4
        for i, ex in enumerate(examples):
            r, c = divmod(i, cols)

            card = ttk.Frame(frame, style="Card.TFrame", padding=10)
            card.grid(row=r, column=c, padx=10, pady=10, sticky="nsew")

            img = load_image(os.path.join(topic_path, ex, "icon.png"))
            if img:
                lbl = ttk.Label(card, image=img)
                lbl.image = img
                lbl.pack()

            ttk.Button(card, text=ex,
                       command=lambda e=ex: self.show_run(topic, e)).pack(fill=tk.X)

        for c in range(cols):
            frame.grid_columnconfigure(c, weight=1)

    # ---------------- Run Screen ---------------- #

    def show_run(self, topic, example):
        self.clear()

        back_cmd = self.show_topics if example is None else lambda: self.show_examples(topic)

        ttk.Button(self.main_frame, text="← Back",
                   command=back_cmd).pack(anchor="w")

        label = topic if example is None else f"{topic} / {example}"

        ttk.Label(self.main_frame, text=label,
                  font=("Arial", 18, "bold")).pack(pady=10)

        self.use_docker = tk.BooleanVar()

        ttk.Checkbutton(self.main_frame, text="Use Docker",
                        variable=self.use_docker).pack(anchor="w")

        # Buttons
        btn_frame = ttk.Frame(self.main_frame)
        btn_frame.pack(pady=10)

        ttk.Button(btn_frame, text="Run Simulation",
                   command=lambda: self.run(topic, example)).pack(side=tk.LEFT, padx=8)

        ttk.Button(btn_frame, text="Open Results",
                   command=lambda: self.open_folder(topic, example)).pack(side=tk.LEFT, padx=8)

        ttk.Button(btn_frame, text="Open ParaView",
                   command=lambda: self.open_paraview(topic, example)).pack(side=tk.LEFT, padx=8)

        # Output
        self.output = scrolledtext.ScrolledText(self.main_frame, height=20)
        self.output.pack(fill=tk.BOTH, expand=True)

    # ---------------- Run Logic ---------------- #

    def run(self, topic, example):
        self.output.delete(1.0, tk.END)
        threading.Thread(target=self._run_thread,
                         args=(topic, example), daemon=True).start()

    def _run_thread(self, topic, example):
        if example is None:
            script_path = topic
        else:
            script_path = f"{topic}/{example}"

        mode = "docker" if self.use_docker.get() else "local"

        try:
            process = subprocess.Popen(
                ["bash", "run.sh", script_path, mode],
                cwd=BASE_DIR,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )

            for line in process.stdout:
                self.safe_print(line)

            process.wait()

            if process.returncode != 0:
                self.safe_print("\n❌ Failed\n")
            else:
                self.safe_print("\n Done\n")

        except Exception as e:
            self.safe_print(f"Error: {e}\n")

    # ---------------- Open Folder ---------------- #

    def open_folder(self, topic, example):
        if example is None:
            base_folder = os.path.join(BASE_DIR, topic)
        else:
            base_folder = os.path.join(BASE_DIR, topic, example)

        # detect results folder
        for name in ["results", "output", "outputs"]:
            candidate = os.path.join(base_folder, name)
            if os.path.isdir(candidate):
                folder = candidate
                break
        else:
            folder = base_folder

        try:
            if "microsoft" in platform.uname().release.lower():
                win_path = subprocess.check_output(
                    ["wslpath", "-w", folder]
                ).decode().strip()
                subprocess.run(["explorer.exe", win_path])

            elif platform.system() == "Windows":
                os.startfile(folder)

            elif platform.system() == "Darwin":
                subprocess.run(["open", folder])

            else:
                subprocess.run(["xdg-open", folder])

        except:
            messagebox.showinfo("Path", folder)

    # ---------------- ParaView ---------------- #

    def open_paraview(self, topic, example):
        if example is None:
            folder = os.path.join(BASE_DIR, topic)
        else:
            folder = os.path.join(BASE_DIR, topic, example)

        pvsm = [f for f in os.listdir(folder) if f.endswith(".pvsm")]

        if not pvsm:
            messagebox.showerror("Error", "No .pvsm file found")
            return

        pvsm_path = os.path.join(folder, pvsm[0])

        try:
            if "microsoft" in platform.uname().release.lower():
                win_path = subprocess.check_output(
                    ["wslpath", "-w", pvsm_path]
                ).decode().strip()

                subprocess.run(["cmd.exe", "/c", "start", "", win_path])

            else:
                subprocess.run(["paraview", pvsm_path])

        except Exception as e:
            messagebox.showerror("Error", str(e))


# ---------------- MAIN ---------------- #

def main():
    if GUI_MODE:
        ExampleSelector(root)
        root.mainloop()


if __name__ == "__main__":
    main()