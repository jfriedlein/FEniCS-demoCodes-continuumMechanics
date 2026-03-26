# powered by perplexity.ai

import os
import tkinter as tk
from tkinter import messagebox
from PIL import Image, ImageTk
import re

# requires pip install Pillow

BASE_DIR = "./"  # script should be called from the main folder (where "run.sh" and "exampleSelector.py" are located)


def load_image(path, max_width=90, max_height=90):
    """Load and resize an image, preserving aspect ratio, or return None on error."""
    try:
        img = Image.open(path)
        img.thumbnail((max_width, max_height), Image.LANCZOS)
        return ImageTk.PhotoImage(img)
    except Exception as e:
        #print(f"Could not load image {path}: {e}")
        return None


class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Topic & Example Selector")
        self.root.geometry("700x500")

        # Ensure the base dir exists
        if not os.path.exists(BASE_DIR):
            messagebox.showerror(
                "Error",
                f"Base directory '{BASE_DIR}' not found.\n"
                "Please create it as described in the comments."
            )
            root.quit()
            return

        # Start with the topic grid
        self.show_topics()

    def clear_window(self):
        """Remove all widgets from the window."""
        for widget in self.root.winfo_children():
            widget.destroy()

    def show_topics(self):
        """Show a grid of topic buttons (only folders starting with '##-')."""
        self.clear_window()
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        # --- Title bar ---
        title_frame = tk.Frame(self.root, pady=10)
        title_frame.grid(row=0, column=0, sticky="nwe")

        title_label = tk.Label(
            title_frame,
            text="Select a topic",
            font=("Helvetica", 14, "bold"),
        )
        title_label.pack(padx=10)

        # --- Topics grid (no scroll) ---
        main_frame = tk.Frame(self.root, pady=10)
        main_frame.grid(row=1, column=0, sticky="nsew")

        # Configure resizing so columns share space
        cols = 3  # fewer columns → wider but no scrollbar
        for c in range(cols):
            main_frame.grid_columnconfigure(c, weight=1)

        # --- Filter only folders that start with two digits + dash ---
        topics_raw = [
            name
            for name in os.listdir(BASE_DIR)
            if os.path.isdir(os.path.join(BASE_DIR, name))
        ]

        pattern = re.compile(r"^\d{2}-")
        topics = [name for name in topics_raw if pattern.match(name)]

        # --- Arrange topic buttons in grid ---
        for idx, topic in enumerate(topics):
            row = idx // cols
            col = idx % cols

            frame = tk.Frame(
                main_frame,
                relief="raised",
                borderwidth=2,
                padx=5,
                pady=5,
            )
            frame.grid(row=row, column=col, padx=5, pady=5, sticky="ew")

            # Button adapts width to text
            btn = tk.Button(
                frame,
                text=topic,
                command=lambda t=topic: self.show_examples(t),
            )
            btn.pack(expand=True, fill="x", pady=5)

    def show_examples(self, topic):
        """Show examples (subfolders) for the given topic."""
        topic_path = os.path.join(BASE_DIR, topic)
        self.clear_window()
        self.root.grid_rowconfigure(0, weight=0)  # title row
        self.root.grid_rowconfigure(1, weight=1)  # examples grid
        self.root.grid_columnconfigure(0, weight=1)

        # --- Title bar: selected topic ---
        title_frame = tk.Frame(self.root, pady=10)
        title_frame.grid(row=0, column=0, sticky="nwe")

        title_label = tk.Label(
            title_frame,
            text=f"Topic: {topic}",
            font=("Helvetica", 14, "bold"),
        )
        title_label.pack(side="left", padx=10)

        # Back button to return to topic selection
        back_btn = tk.Button(
            title_frame,
            text="← Back",
            command=self.show_topics,
        )
        back_btn.pack(side="right", padx=10)

        # --- Examples grid (no scroll) ---
        main_frame = tk.Frame(self.root, pady=10)
        main_frame.grid(row=1, column=0, sticky="nsew")

        cols = 3  # fewer columns, better fit in window
        for c in range(cols):
            main_frame.grid_columnconfigure(c, weight=1)

        examples = [
            name
            for name in os.listdir(topic_path)
            if os.path.isdir(os.path.join(topic_path, name))
        ]

        for idx, example in enumerate(examples):
            row = idx // cols
            col = idx % cols

            frame = tk.Frame(
                main_frame,
                relief="raised",
                borderwidth=2,
                padx=5,
                pady=5,
            )
            frame.grid(row=row, column=col, padx=5, pady=5, sticky="ew")

            icon_path = os.path.join(topic_path, example, "icon.png")
            img = load_image(icon_path, max_width=90, max_height=90)
            if img:
                img_label = tk.Label(frame, image=img)
                img_label.image = img  # Keep a reference
                img_label.pack(pady=2)

            btn = tk.Button(
                frame,
                text=example,
                command=lambda ex=example, t=topic: self.select_example(t, ex),
            )
            btn.pack(expand=True, fill="x", pady=2)

    def select_example(self, topic, example):
        """Callback when an example is clicked; outputs the relative path."""
        relative_path = os.path.join(topic, example).replace("\\", "/")  # normalize
        #messagebox.showinfo(
        #    "Selected Path",
        #    f"Selected relative path:\n{relative_path}"
        #)
        #print(f"Selected relative path: {relative_path}")
        print(relative_path)      # <-- this is what Bash will capture
        self.root.quit()          # close GUI


def main():
    root = tk.Tk()
    app = App(root)
    root.mainloop()


if __name__ == "__main__":
    main()

