import tkinter as tk
import tkinter.ttk as ttk
import sys


BG = "#edecd6"
FG = "grey15"


class SpriteEditor(tk.Tk):
    def __init__(self, output, size, px):
        super().__init__()
        self.output = output
        self.size = size
        self.px = px
        self.ncells = size // px
        self.titletxt = f"{output} [{self.ncells}x{self.ncells}]"
        self.title(self.titletxt)

        self.resizable(False, False)

        self.frame = ttk.Frame(self)
        self.frame.pack(fill="both", expand=True, padx=10, pady=10)

        self.canvas = tk.Canvas(
            self.frame,
            width=size,
            height=size,
            bg=BG,
            borderwidth=0,
            highlightthickness=0,
        )
        self.canvas.pack(fill="both", expand=True)

        self.canvas.bind("<B1-Motion>", self.paint_fg)
        self.canvas.bind("<B3-Motion>", self.paint_bg)
        self.canvas.bind("<Motion>", self.show_position)

        self.buttons = ttk.Frame(self.frame)
        self.buttons.grid_columnconfigure(0, weight=1)
        self.buttons.grid_columnconfigure(1, weight=1)
        self.buttons.pack(fill="both", expand=True)
        self.clear_button = ttk.Button(
            self.buttons, text="Borrar todo", command=self.clear_canvas
        )
        self.clear_button.grid(row=0, column=0, padx=(0, 5), pady=(10, 0), sticky="ew")
        self.save_button = ttk.Button(
            self.buttons, text="Guardar", command=self.save_ppm
        )
        self.save_button.grid(row=0, column=1, padx=(5, 0), pady=(10, 0), sticky="ew")

        self.pixels = [
            [0 for _ in range(self.ncells)]
            for _ in range(self.ncells)
        ]

    def is_inside_canvas(self, x, y):
        return 0 <= x < self.size and 0 <= y < self.size

    def paint_fg(self, event):
        x, y = event.x // self.px, event.y // self.px
        if self.is_inside_canvas(event.x, event.y):
            self.pixels[y][x] = 1
            self.draw_pixel(x, y, FG)
            self.show_position(event)

    def paint_bg(self, event):
        x, y = event.x // self.px, event.y // self.px
        if self.is_inside_canvas(event.x, event.y):
            self.pixels[y][x] = 0
            self.draw_pixel(x, y, BG)
            self.show_position(event)

    def draw_pixel(self, x, y, color):
        x1, y1 = x * self.px, y * self.px
        x2, y2 = x1 + self.px, y1 + self.px
        self.canvas.create_rectangle(x1, y1, x2, y2, fill=color, outline="")

    def show_position(self, event):
        x, y = event.x // self.px, event.y // self.px
        self.title(f"{self.titletxt} ({x}, {y})")

    def clear_canvas(self):
        self.canvas.delete("all")
        self.canvas.config(bg=BG)
        self.pixels = [
            [0 for _ in range(self.ncells)]
            for _ in range(self.ncells)
        ]

    def save_ppm(self):
        with open(f"{OUTPUT}.pbm", "w") as file:
            file.write("P1\n")
            file.write(f"{self.ncells} {self.ncells}\n")
            for row in self.pixels:
                file.write(f"{' '.join([str(p) for p in row])}\n")
        print(f"Perfil guardado en {OUTPUT}.pbm")


def _show_help():
    print("Error parámetros, uso correcto:")
    print("    profile_designer.exe nombre_perfil numero_celdas [tamaño_lienzo]")
    print("\nEjemplo:")
    print(f"    profile_designer.exe mi_perfil 10")
    sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        OUTPUT = sys.argv[1]
    else:
        _show_help()
    if len(sys.argv) > 2 and sys.argv[2].isdigit():
        CELLS = int(sys.argv[2])
    else:
        _show_help()

    if len(sys.argv) > 3:
        if sys.argv[3].isdigit():
            SIZE = int(sys.argv[3])
        else:
            _show_help()
    else:
        SIZE = 400

    PX = round(SIZE / CELLS)
    SIZE = PX * CELLS

    app = SpriteEditor(OUTPUT, SIZE, PX)
    app.mainloop()
