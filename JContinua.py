# -*- coding: utf-8 -*-
#
#  JContinua - programa para ejercitar el trazado de diagramas de esfuerzos 
#  y desplazamientos en vigas.
#  
#  Copyright 2021: Juan Carlos del Caño Sánchez, profesor en la Escuela de
#  Ingenierías Industriales de Valladolid, España (actualmente retirado).
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

from tkinter import *
from tkinter import ttk, messagebox, filedialog

import numpy as np
import numpy.polynomial.polynomial as poly

import os
from os import path, kill, getppid, _exit

from random import randrange as rr
from copy import deepcopy
import signal


class Nodo:
    # tiene un cuadro(=tk.Frame) para hacer grid, tiene 3 botones 
    # que se autocomprueban cuando cambias uno, tiene lista de 3 BooleanVar
    # (bv_gdl) para los 3 botones, tiene asociado un objeto union que se 
    # asignara al pulsar "mostrar" y que tendra el dibujo oportuno

    def verif_shareux(self):
        if not self.bv_gdl[0].get() and not self.bv_gdl[1].get(): 
            self.bv_gdl[0].set(True)

    def verif_shareuy(self):
        if not self.bv_gdl[0].get() and not self.bv_gdl[1].get(): 
            self.bv_gdl[1].set(True)

    def __init__(self, inodo, x=0.):
        
        self.x=x
        self.cuadro= ttk.Frame(frame_general)
        
        self.bv_gdl=[BooleanVar(value=True), 
                     BooleanVar(value=True), 
                     BooleanVar(value=True)]

        self.lbln= ttk.Label(self.cuadro, text=str(inodo), background='white')
        
        self.btn_shareux= ttk.Checkbutton(self.cuadro, style='jc_blue.TCheckbutton',                 variable=self.bv_gdl[0], command=self.verif_shareux)
        self.btn_shareuy= ttk.Checkbutton(self.cuadro, style='jc_blue.TCheckbutton',    
                         variable=self.bv_gdl[1], command=self.verif_shareuy)
        self.btn_sharefi= ttk.Checkbutton(self.cuadro, style='jc_blue.TCheckbutton',
                         variable=self.bv_gdl[2] )

        self.btn_shareux.state(['!alternate','selected'])
        self.btn_shareuy.state(['!alternate','selected'])
        self.btn_sharefi.state(['!alternate','selected'])
        
        self.lbln.grid(row=0, column=0)
        self.btn_shareux.grid(row=1, column=0)
        self.btn_shareuy.grid(row=2, column=0)
        self.btn_sharefi.grid(row=3, column=0)





class Tramo:
    def __init__(self, ctes, inodos):
        self.inodos= inodos # lista con sus dos codigos de nodo
        self.ctes= ctes # lista [L,E,I,A]
        self.seti=Setde_s() # condiciones de apoyo en extremos
        self.setj=Setde_s()
        self.carga=[[0., 0.], [0., 0.], [0., 0.], [0., 0.], [0., 0.]] # px py fx fy mz 
        self.fa_tr = [ [0,0,0,0,0,0], [0,0,0,0,0,0] ] # V,M,fi,u,N,w, de ini & de fin (fuerzas!)
        self.funciones=[] # funciones(x) de idem V,M,fi,u,N,w
        self.haycero=[ [],[],[],[],[],[] ] # vacios ó ceros interiores al tramo de esas funciones
        self.maximos=[ [],[],[],[],[],[] ] # vacios o con duplas [x_max,valor]
        self.reacc=[  [ [],[],[] ]  ,  [ [],[],[] ]  ] # iój / V,MóN / vacio ó valor 

    def ponfunciones(self):
        V_x = poly.Polynomial([ -self.fa_tr[0][0], - self.carga[1][0], 
                - (self.carga[1][1]-self.carga[1][0])/(2*self.ctes[0])  ])
        M_x = V_x.integ()
        M_x.coef[0]= -self.fa_tr[0][1]
        fi_x = M_x.integ() / (self.ctes[1]*self.ctes[2]) 
        fi_x.coef[0]=self.fa_tr[0][2]
        u_x = -fi_x.integ() 
        u_x.coef[0] = self.fa_tr[0][3]
        
        r_x = poly.Polynomial([ self.carga[0][0], (self.carga[0][1]-self.carga[0][0])/self.ctes[0]])
        N_x= -r_x.integ()
        N_x.coef[0]= -self.fa_tr[0][4]
        w_x= N_x.integ()/ (self.ctes[1]*self.ctes[3])
        w_x.coef[0]= self.fa_tr[0][5]
        
        self.funciones=[V_x, M_x, fi_x, u_x, N_x, w_x]


    def pinta_py(self, lienzo, ypx0, scal_py, acota):
        f= lambda x: ( self.carga[1][0] + 
                    x*(self.carga[1][1]-self.carga[1][0])/self.ctes[0] )
        i, kk = 0, [] # para linea superior
        for x in np.linspace(0, self.ctes[0], 4):
            y=f(x)
            x_pixel = pixx(nodos[self.inodos[0]].x + x)
            y_pixel_ini=ypx0- round( y*scal_py )
            y_pixel_fin=ypx0
            lienzo.create_line(x_pixel, y_pixel_ini, x_pixel, y_pixel_fin, 
                arrow='last', width=2, fill='#FF3E3E', arrowshape=(7,11,4),
                tags=('carga'))
            if i in (0,3): kk.append([x_pixel,y_pixel_ini])
            i += 1
        
        lienzo.create_line(kk[0][0], kk[0][1], kk[1][0], kk[1][1],
                width=1, fill='#FF3E3E', tags=('carga'))
        if acota:
            c='sw' if self.carga[1][0]>0 else 'nw'
            lienzo.create_text(kk[0][0], kk[0][1], text='{:.3g}'.format(abs(self.carga[1][0])), 
                                fill='red', anchor=c)
            c='se' if self.carga[1][1]>0 else 'ne'
            lienzo.create_text(kk[1][0]+rr(-7,-3), kk[1][1]+rr(3,7), 
                        text='{:.3g}'.format(abs(self.carga[1][1])), fill='red', anchor=c)


    def pinta_px(self, lienzo, ypx0, scal_px, yoffset, acota):
        x0, x1 = pixx(nodos[self.inodos[0]].x)+4, pixx(nodos[self.inodos[1]].x)-4
        ybase= ypx0-yoffset
        y0 = ybase +round(self.carga[0][0]*scal_px)
        y1 = ybase +round(self.carga[0][1]*scal_px) 
        lista_pa_trazar=[x0,ybase, x1,ybase, x1,y1, x0,y0]
        lienzo.create_polygon(lista_pa_trazar, outline='blue', fill='#C1D6EB',
                         stipple='gray12', tags=('carga'))
        if acota:
            c='nw' if self.carga[0][0]>0 else 'sw'
            lienzo.create_text(x0,round(y0), text='{:.3g}'.format(abs(self.carga[0][0])), 
                                fill='blue', anchor=c)
            c='ne' if self.carga[0][1]>0 else 'se'
            lienzo.create_text(x1,round(y1), text='{:.3g}'.format(abs(self.carga[0][1])), 
                                fill='blue', anchor=c)

        # para unas flechicas:
        signoi= np.sign(self.carga[0][0]) or np.sign(self.carga[0][1]) # si es 0 se pone el del otro
        signoj= np.sign(self.carga[0][1]) or np.sign(self.carga[0][0]) # idem. 
                                                        # Pueden quedar ambos =0 pero no uno solo
        a= round((x1-x0)/8)
        if signoi:
            x0, x1= x0+a, x0+2*a
            if signoi==-1: x0,x1 = x1,x0
            b=15 if signoi>0 else -15
            y0= ybase + b
            lienzo.create_line(x0,y0,x1,y0, arrow='last',width=2, fill='blue',
                        arrowshape=(7,11,4), tags=('carga'))
        # if signoj: no hace falta, si signoi!=0 -> signoj tb lo sera
            x1 = pixx(nodos[self.inodos[1]].x)-25
            x0 = x1-a
            if signoj==-1: x0,x1 = x1,x0
            b=15 if signoj>0 else -15
            y0= ybase + b
            lienzo.create_line(x0,y0,x1,y0, arrow='last',width=2, fill='blue',
                        arrowshape=(7,11,4), tags=('carga'))


    def pinta_fy(self, lienzo, ypx0, scal_fy, acota):
        
        if self.carga[3][0]: # F en nodo izdo
            nodo= nodos[self.inodos[0]]
            x_pixel=pixx(nodo.x)
            if not nodo.bv_gdl[1].get(): x_pixel += 12
            y_pixel_ini= ypx0-20- abs(round(self.carga[3][0]*scal_fy))
            y_pixel_fin= ypx0-20
            if self.carga[3][0] < 0.0:
                y_pixel_ini, y_pixel_fin = y_pixel_fin, y_pixel_ini
            lienzo.create_line(x_pixel, y_pixel_ini, x_pixel, y_pixel_fin, 
                    arrow='last', width=6, fill='#FFB4A6', arrowshape=(16,26,12),
                    tags=('carga'))
            if acota:
                y= round(y_pixel_ini+ y_pixel_fin)/2
                lienzo.create_text(x_pixel, y,  font=('Times','12','bold'), fill='#FFB4A6',
                                    text='{:.3g}'.format(abs(self.carga[3][0])), anchor='w')

        if self.carga[3][1]: # F en nodo dcho
            nodo= nodos[self.inodos[1]]
            x_pixel=pixx(nodo.x)
            if not nodo.bv_gdl[1].get(): x_pixel -= 12
            y_pixel_ini= ypx0-20- abs(round(self.carga[3][1]*scal_fy))
            y_pixel_fin= ypx0-20
            if self.carga[3][1] < 0.0:
                y_pixel_ini, y_pixel_fin = y_pixel_fin, y_pixel_ini
            lienzo.create_line(x_pixel, y_pixel_ini, x_pixel, y_pixel_fin, 
                    arrow='last', width=6, fill='#FFB4A6', arrowshape=(16,26,12),
                    tags=('carga'))
            if acota:
                y= round(y_pixel_ini+ y_pixel_fin)/2
                lienzo.create_text(x_pixel, y, font=('Times','12','bold'), fill='#FFB4A6',
                                    text='{:.3g}'.format(abs(self.carga[3][1])), anchor='e')


    def pinta_fx(self, lienzo, ypx0, scal_fx, acota):
        
        if self.carga[2][0]: # F en nodo izdo
            nodo= nodos[self.inodos[0]]
            c= round(self.carga[2][0]*scal_fx)
            if np.isclose(nodo.x, 0.) and c<0:
                x_pixel_ini = pixx(nodo.x) + abs(c)
                x_pixel_fin = pixx(nodo.x)
            else:
                x_pixel_ini, x_pixel_fin = pixx(nodo.x), pixx(nodo.x)+c
            if not nodo.bv_gdl[0].get(): 
                x_pixel_ini += 12
                x_pixel_fin += 12
            
            y_pixel = ypx0 - 20
            lienzo.create_line(x_pixel_ini, y_pixel, x_pixel_fin, y_pixel, 
                    arrow='last', width=4, fill='#98BEE2', arrowshape=(14,22,10),
                    tags=('carga'))
            if acota:
                c= (x_pixel_fin+ x_pixel_ini )/2
                lienzo.create_text(c, y_pixel+15,font=('Times','12','bold'), fill='#98BEE2', text='{:.3g}'.format(abs(self.carga[2][0])), anchor='center')

        if self.carga[2][1]: # F en nodo dcho
            nodo= nodos[self.inodos[1]]
            x_pixel_ini=pixx(nodo.x)
            if not nodo.bv_gdl[0].get(): x_pixel_ini -= 12
            x_pixel_fin= x_pixel_ini +round(self.carga[2][1]*scal_fx)
            y_pixel = ypx0 - 40
            lienzo.create_line(x_pixel_ini, y_pixel, x_pixel_fin, y_pixel, 
                    arrow='last', width=4, fill='#98BEE2', arrowshape=(14,22,10),
                    tags=('carga'))
            if acota:
                c= (x_pixel_fin+ x_pixel_ini )/2
                lienzo.create_text(c, y_pixel+15,font=('Times','12','bold'), fill='#98BEE2', text='{:.3g}'.format(abs(self.carga[2][1])), anchor='center')


    def pinta_mz(self, lienzo, ypx0, acota):
        
        def momentito(x_pixel, xoffset, valor, acota, ladocota):
            # dibuja un momento 0-180º en x_pixel+-25, y=290-+25, con flecha
            # a la izda si signo=1. El xoffset será 0 ó +-12 
            signo= np.sign(valor)
            a, b = [x_pixel-25+xoffset, ypx0-35], [x_pixel+25+xoffset, ypx0+15]
            lienzo.create_arc(a[0],a[1],b[0],b[1], outline='orange',
                        style='arc', start=0, extent=180, width=3,tags=('carga'))
            c= a[0]+1 if signo >0 else b[0]-1
            lienzo.create_line(c, ypx0-8, c, ypx0-5, 
                    arrow='last', width=3, fill='orange', arrowshape=(8,13,6),
                    tags=('carga'))
            if acota:
                d=40 if ladocota=='dcha' else -40
                lienzo.create_text(x_pixel+xoffset+d, ypx0-15,font=('Times','12','bold'), fill='orange', text='{:.3g}'.format(abs(valor)))
        
        if self.carga[4][0]: # M en nodo izdo
            nodo= nodos[self.inodos[0]]
            x_pixel=pixx(nodo.x)
            if nodo.bv_gdl[2].get(): # gdl compartido, dibujo arriba
                momentito(x_pixel, 0, self.carga[4][0], acota, 'dcha')
            else: # gdl no compartido, dibujo a la dcha
                momentito(x_pixel, 12, self.carga[4][0], acota, 'dcha')

        if self.carga[4][1]: # M en nodo dcho
            nodo= nodos[self.inodos[1]]
            x_pixel=pixx(nodo.x)
            if nodo.bv_gdl[2].get(): # gdl compartido, dibujo arriba lo mismo
                momentito(x_pixel, 0, self.carga[4][1], acota, 'izda')
            else: # gdl no compartido, dibujo a la dcha
                momentito(x_pixel, -12, self.carga[4][1], acota, 'izda')



class Setde_s:
    # la s final es de sustentacion "set-de-s". 
    # Hare grid del cuadro. Los valores estan en self.bv_restri (BoolVar), no tengo
    # que manipular el valor del bptón directamente con cosas como:
    # btn_ux.instate(['selected']) -> True o False 
    # btn_ux.state() ->  una tupla tipo (['selected'])

    def __init__(self):
        self.cuadro= ttk.Frame(frame_general,style='jc_green.TFrame')
        
        self.bv_restri=[BooleanVar(value=False), BooleanVar(value=False), 
                        BooleanVar(value=False), BooleanVar(value=False), 
                        BooleanVar(value=False)]
        
        self.btn_ux= ttk.Checkbutton(self.cuadro, style='jc_green.TCheckbutton',
                            variable=self.bv_restri[0])
        self.btn_uy= ttk.Checkbutton(self.cuadro, style='jc_green.TCheckbutton',
                            variable=self.bv_restri[1])
        self.btn_fi= ttk.Checkbutton(self.cuadro, style='jc_green.TCheckbutton', 
                            variable=self.bv_restri[2])
        self.btn_ky= ttk.Checkbutton(self.cuadro, style='jc_green.TCheckbutton',
                            variable=self.bv_restri[3],)
        self.btn_kt= ttk.Checkbutton(self.cuadro, style='jc_green.TCheckbutton',
                            variable=self.bv_restri[4])

        self.btn_ux.grid(row=0, column=0, sticky='e')
        self.btn_uy.grid(row=1, column=0, sticky='e')
        self.btn_fi.grid(row=2, column=0, sticky='e')
        self.btn_ky.grid(row=3, column=0, sticky='e')
        self.btn_kt.grid(row=4, column=0, sticky='e')
        
        self.valores=[0.0, 0.0, 0.0, 0.0, 0.0]
        self.entries=[ ttk.Entry(frame_general, width=3),
                ttk.Entry(frame_general, width=3),ttk.Entry(frame_general, width=3),
                ttk.Entry(frame_general, width=3),ttk.Entry(frame_general, width=3)]


    def grid_un_entry_s(self, i=0, xnodo=0, izda=False):
        # lo hago en el canvas (no en cuadro) para que tenga tag
        self.entries[i].delete(0,'end')
        self.entries[i].insert(0,str(self.valores[i]))
        
        px= pixx(xnodo)-16 if izda else pixx(xnodo)+16
        py= 400+ i*21
        frame_lienzo.create_window(px,py, anchor='n', window=self.entries[i], 
                tags=('entry_s'))
    
    
    def desempareja(self):
        self.bv_restri=[BooleanVar(value=False), BooleanVar(value=False), 
                        BooleanVar(value=False), BooleanVar(value=False), 
                        BooleanVar(value=False)]
        self.btn_ux.configure(variable=self.bv_restri[0])
        self.btn_uy.configure(variable=self.bv_restri[1])
        self.btn_fi.configure(variable=self.bv_restri[2])
        self.btn_ky.configure(variable=self.bv_restri[3])
        self.btn_kt.configure(variable=self.bv_restri[4])



class ToolTip(object):
    # para los cuadradillos informativos emergentes
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()



def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)



def ayuda():

    v5=Toplevel(v0)
        
    texto='''- Licencia y presentación.
    
JContinua es un programa ofrecido bajo licencia GPL de Gnu. Básicamente usted puede estudiarlo, usarlo y distribuirlo libremente, pero si realiza modificaciones al programa y quiere hacerlas públicas, debe hacerlo bajo la misma licencia GPL y citando al autor original. Puede consultar los términos completos de dicha licencia en http://www.gnu.org/licenses/gpl.html

El programa realiza el trazado de los diagramas de esfuerzos y desplazamientos en vigas, tal como se suelen estudiar en la mayoría de titulaciones de arquitectura e ingeniería. Tiene el formato de un programa de cálculo pero su finalidad principal es didáctica: se espera que el estudiante lo use como una herramienta de apoyo y comprobación en sus ejercicios de trazado de diagramas. Las categorías de conceptos (uniones, apoyos, cargas) se presentan claramente estructuradas, lo que favorecerá una correcta comprensión del modelo de cálculo.

Esta es la primera versión de JContinua que se hace pública, y su estado de desarrollo es considerado "funcional sin extras". La conveniencia de hacerlo disponible en una fecha concreta ha motivado que se aplace la incorporación de algunas características no esenciales.

- Uso del programa.

El programa asume que los datos están suministrados coherentemente en un cierto sistema de unidades elegido por el usuario. Los resultados se expresan en esas mismas unidades.

Se comienza pulsando el botón "tramos" en donde se especifican las longitudes de cada tramo (hasta 8 tramos, deje en blanco los no usados), su módulo de Young (E), momento de inercia (I) y el área de la sección. Use un número de tramos adecuado teniendo en cuenta que cualquier discontinuidad ha de ocurrir en un nodo o extremo de tramo. Es el caso de una fuerza o momento concentrados, un apoyo, el comienzo o el fin de una carga distribuida, etc.

Seguidamente se definen las uniones entre tramos (en azul). Por defecto aparecen compartidos tanto el movimiento horizontal (ux) como el vertical (uy) como el giro de la sección (fi), lo que corresponde a continuidad de la viga en esa unión. Puede incorporar libertades desmarcando los botones azules correspondientes. No se permite desmarcar simultáneamente ux & uy por resultar un tipo de unión difícilmente materializable. Se pulsa el botón "incorporar uniones" al terminar esta etapa, con lo que se dibujan los símbolos correspondientes en la viga.

Seguidamente se definen los apoyos y resortes que forman la sustentación (botones en verde). Se admiten apoyos que impongan ux, uy, fi, asi como resortes de desplazamiento vertical (ky) y de giro (kt). Una vez marcados los botones requeridos aparecen los símbolos de sustentación correspondientes al pulsar "incorporar apoyos". 

En el caso en que algún apoyo tenga desplazamiento ó giro conocido ("asiento"), o que existan resortes, necesitaremos pulsar el botón "valores no nulos" para especificar el valor del asiento o el valor de la constante de resorte. Por favor lea la nota emergente que aparece al pasar el ratón sobre este botón. Si cambia algún apoyo posteriormente, los valores no nulos se pondrán de nuevo a cero.

Finalmente se definen las cargas pulsando el botón "definir cargas". Se admiten cargas perpendiculares a la viga tanto concentradas (Fy) como distribuidas (py), asi como momentos concentrados perpendiculares al plano del dibujo (Mz). También se admiten cargas colineales con la barra, tanto concentradas (Fx) como distribuidas (px). Las cargas distribuidas pueden ser constantes o lineales en cada tramo. Por favor lea la nota que emerge al pasar el ratón sobre el cuerpo de la ventana, que indica cómo especificar las cargas concentradas cuando el movimiento correspondiente esté (o no) compartido en la unión. Al pulsar el botón "hecho" en esta ventana la misma se cierra y se dibujan las cargas en colores cálidos (cercanos al rojo).

En este momento tenemos el modelo de la viga completamente especificado. Podemos pulsar el botón "calcular". Se generará un trazado con los diagramas de esfuerzos cortantes, momentos flectores, giros y desplazamientos, y otro trazado con los diagramas de esfuerzo axil y desplazamiento axil. También se representan las reacciones de la sustentación (en color magenta). En la salida de texto se especifican con más decimales los valores acotados en las gráficas.

El botón "imprimir" genera dos ficheros postcript y un fichero de texto en el directorio deseado, correspondiendo a cada una de las ventanas gráficas obtenidas y a la salida de texto. 
Nota: debido a un error en el backend utilizado, se requiere que la cabecera de las ventanas gráficas esté visible en pantalla en el momento de "imprimir" para que se guarden completas.

- Convenios de signos.

Las magnitudes vectoriales contenidas en el plano del dibujo son positivas hacia abajo ó hacia la derecha. Es el caso de Fy, Fx, py, px, uy, ux. También de las reacciones de fuerza, y las fuerzas en los extremos de tramo que proporciona la salida en texto de resultados. 

Las magnitudes vectoriales perpendiculares al plano del dibujo son positivas en sentido saliente del plano del mismo. Es el caso de los momentos concentrados aplicados, las reacciones de momento, el giro fi de la sección, y los momentos en los extremos de tramo que proporciona la salida en texto de resultados. 

Nótese que en la salida gráfica, tanto las magnitudes vectoriales (reacciones...) como las no vectoriales (esfuerzo cortante, momento flector, esfuerzo axil) se representan en  en la salida gráfica "como son" y con cotas sin signo tal como indica la buena práctica. En particular, los esfuerzos no aparecen en la salida de texto por lo que el usuario no necesita saber el convenio de signos con el que se han calculado. Como curiosidad, los diagramas que se trazan en la parte inferior de las gráficas corresponden a esfuerzos, giros ó desplazamientos que internamente se consideran positivos, siendo negativos los que se trazan en la parte superior.

            ¡Espero que JContinua le sea útil!
             Juan Carlos del Caño
             Profesor en la EII de Valladolid (retirado)

            ____________________________________
'''


    ttk.Label(v5, text='Notas breves',
        font=('', 16)).grid(row=0, column=0, columnspan=2, pady=5)

    tcaja = Text(v5, width=45, height=30,wrap='word', font=('Sans',9),
        background='#D9D9D9', foreground='green', border=None, padx=20, pady=12)
    tcaja.grid(column=0, row=1, padx=8, sticky=(N,W,E,S))
    tcaja.insert('1.0',texto)

    scb = ttk.Scrollbar(v5, orient=VERTICAL, command=tcaja.yview)
    scb.grid(column=1, row=1, sticky='ns')
    
    tcaja['yscrollcommand'] = scb.set

    ttk.Button(v5, text='Entendido', width=9, command=v5.destroy).grid(
        column=0, row=2, pady=4, columnspan=2)

    tcaja['state']='disabled'

    v5.grid_columnconfigure(0, weight=1)
    v5.grid_rowconfigure(0, weight=1)
    v5.grid_rowconfigure(1, weight=4)
    v5.grid_rowconfigure(2, weight=1)
    v5.geometry('+250+70')

    v5.focus()
    v5.mainloop()



def salir():
    global v3,v4
    if messagebox.askokcancel(message='¿Quiere cerrar el programa? ',default='cancel',
            icon='question', title='Confirmacion:',parent=v0) :
        cerrar_graf()
        try:
            v5.destroy()
        except:
            pass
        v0.destroy()
        try: # si es linux que cierre el terminal tambien
            kill(getppid(), signal.SIGHUP)        
        except:
            exit()
        # con extension .pyw no saca terminal.


def do_uniones():
    global btn_do_cargas, btn_do_apoyos, tramos, nodos
    frame_lienzo.delete('union')

    pinta_uniones(frame_lienzo, 300)
    
    # resetear apoyos y sus botones, y emparejar los de gdl compartido
    frame_lienzo.delete('apoyo')
    for tramo in tramos:
        
        # habilito los botones:
        tramo.seti.btn_ux.state(['!disabled'])
        tramo.seti.btn_uy.state(['!disabled'])
        tramo.seti.btn_fi.state(['!disabled'])
        tramo.seti.btn_ky.state(['!disabled'])
        tramo.seti.btn_kt.state(['!disabled'])
        
        tramo.setj.btn_ux.state(['!disabled'])
        tramo.setj.btn_uy.state(['!disabled'])
        tramo.setj.btn_fi.state(['!disabled'])
        tramo.setj.btn_ky.state(['!disabled'])
        tramo.setj.btn_kt.state(['!disabled'])
        
        # desemparejo:
        tramo.seti.desempareja()
        tramo.setj.desempareja()
        
        # emparejo:
    for inodo in range(1,nnodos-1):
        nodo = nodos[inodo]
        tramo_izdo, tramo_dcho = tramos[inodo-1], tramos[inodo]  

        if nodo.bv_gdl[0].get():
            tramo_izdo.setj.bv_restri[0]=tramo_dcho.seti.bv_restri[0]
            tramo_izdo.setj.btn_ux.config(variable=tramo_izdo.setj.bv_restri[0])
        if nodo.bv_gdl[1].get():
            tramo_izdo.setj.bv_restri[1]=tramo_dcho.seti.bv_restri[1]
            tramo_izdo.setj.bv_restri[3]=tramo_dcho.seti.bv_restri[3]
            tramo_izdo.setj.btn_uy.config(variable=tramo_izdo.setj.bv_restri[1])
            tramo_izdo.setj.btn_ky.config(variable=tramo_izdo.setj.bv_restri[3])
        if nodo.bv_gdl[2].get():
            tramo_izdo.setj.bv_restri[2]=tramo_dcho.seti.bv_restri[2]
            tramo_izdo.setj.bv_restri[4]=tramo_dcho.seti.bv_restri[4]
            tramo_izdo.setj.btn_fi.config(variable=tramo_izdo.setj.bv_restri[2])
            tramo_izdo.setj.btn_kt.config(variable=tramo_izdo.setj.bv_restri[4])


    btn_do_apoyos['state'] = 'normal'


def pinta_uniones(lienzo, ypx0):
    for nodo in nodos[1:-1]:
        # pintar la union:
        name='n'
        if nodo.bv_gdl[0].get(): name += 'ux'
        if nodo.bv_gdl[1].get(): name += 'uy'
        if nodo.bv_gdl[2].get(): name += 'fi'
        lienzo.create_image(pixx(nodo.x), ypx0, image=figs[name],
                                    anchor='center', tags=('union'))
    


def do_apoyos():
    global nnodos, ntramos, nodos, tramos, scal_px
    global btn_do_nonulos
    
    frame_lienzo.delete('apoyo')
    
    # resetear los nonulos a 0
    for tramo in tramos:
        for i in range(5):
            tramo.seti.valores[i]=0.0
            tramo.setj.valores[i]=0.0

    # comprobacion previa de compatibilidad:
    for tramo in tramos:
        if (  (tramo.seti.bv_restri[1].get() and tramo.seti.bv_restri[3].get()) or
              (tramo.seti.bv_restri[2].get() and tramo.seti.bv_restri[4].get()) or
              (tramo.setj.bv_restri[1].get() and tramo.setj.bv_restri[3].get()) or
              (tramo.setj.bv_restri[2].get() and tramo.setj.bv_restri[4].get()) ):
            texto1 ='Incoherencia en la sustentación.'
            texto2 ='Existe un resorte cuyo grado de libertad está restringido.'
            messagebox.showinfo(parent=v0, title='MAL', message=texto1, detail=texto2)
            return
    
    # la info necesaria esta en: 
    ### nodo.bv_gdl[3]                  -> los gdl compartidos
    ### tramo.seti.bv_restri[5] (&setj) -> gdl apoyados + ky, kt
    # con eso puedo poner los simblos. Los valores numericos los pongo aparte (boton)

    pinta_apoyos(frame_lienzo, 300)

    btn_do_nonulos['state'] = 'normal'
    btn_do_cargas['state'] = 'normal'


def pinta_apoyos(lienzo, ypx0):
    
    textos=['ux','uy','fi', 'ky','kt']
    
    # los que son centradoy=True tienen "i/d" final y NO tendran offsetx
    # los que son centradoy=False tendran anchor='n', y pueden tener offsetx:
    centradoy={'skt':False, 'sky':False, 'sux':True, 'suxfi':True, 'sfi':False,
               'suxuy':False, 'suxuyfi':True, 'suy':False, 'suyfi':False}

    for itramo in range(ntramos):
        tramo=tramos[itramo]
        nodoi, nodoj = nodos[itramo], nodos[itramo+1]
        
        # extremo i:
        
        centradox=True # miro si todos los gdl restringidos estan compartidos:
        for i in range(3):
            if tramo.seti.bv_restri[i].get() and not nodoi.bv_gdl[i].get():
                centradox=False
                break
        if itramo==0: centradox=True # es extremo izdo

        # elijo la figura y la pongo (si hay que poner):
        if ( tramo.seti.bv_restri[0].get() or tramo.seti.bv_restri[1].get() or
             tramo.seti.bv_restri[2].get() ): 
            nom_fich='s'
            for i in range(3): # para ux,uy,fi
                if tramo.seti.bv_restri[i].get():
                    nom_fich += textos[i]
            if centradoy[nom_fich]:
                nom_fich += 'i'
                lienzo.create_image(pixx(nodoi.x), ypx0, image=figs[nom_fich], 
                                        tags=('apoyo'))
            else:
                offsetx=0 if centradox else 12
                lienzo.create_image(pixx(nodoi.x)+offsetx, ypx0,
                        image=figs[nom_fich], anchor='n', tags=('apoyo'))

        # si hay que poner resorte (quiza adicional) lo hago:
        if tramo.seti.bv_restri[3].get():
            centradox=True if (nodoi.bv_gdl[1].get() or itramo==0) else False
            offsetx=0 if centradox else 12
            nom_fich='sky'
            lienzo.create_image(pixx(nodoi.x)+offsetx, ypx0, 
                        image=figs[nom_fich], anchor='n', tags=('apoyo'))
        if tramo.seti.bv_restri[4].get():
            centradox=True if (nodoi.bv_gdl[2].get() or itramo==0) else False
            offsetx=0 if centradox else 12
            nom_fich='skt'
            lienzo.create_image(pixx(nodoi.x)+offsetx, ypx0, 
                        image=figs[nom_fich], anchor='n', tags=('apoyo'))

        # extremo j:
        
        centradox=True # miro si todos los gdl restringidos estan compartidos:
        for i in range(3):
            if tramo.setj.bv_restri[i].get() and not nodoj.bv_gdl[i].get():
                centradox=False
                break
        if itramo == ntramos-1: centradox=True # es extremo dcho

        # elijo la figura y la pongo (si hay que poner):
        if ( tramo.setj.bv_restri[0].get() or tramo.setj.bv_restri[1].get() or
             tramo.setj.bv_restri[2].get() ) :
            nom_fich='s'
            for i in range(3): # para ux,uy,fi
                if tramo.setj.bv_restri[i].get():
                    nom_fich += textos[i]
            if centradoy[nom_fich]:
                nom_fich += 'd'
                lienzo.create_image(pixx(nodoj.x), ypx0, image=figs[nom_fich], 
                                        tags=('apoyo'))
            else:
                offsetx=0 if centradox else -12
                lienzo.create_image(pixx(nodoj.x)+offsetx, ypx0,
                        image=figs[nom_fich], anchor='n', tags=('apoyo'))

        # si hay que poner resorte (quiza adicional) lo hago:
        if tramo.setj.bv_restri[3].get():
            centradox=True if (nodoj.bv_gdl[1].get() or itramo== ntramos-1) else False
            offsetx=0 if centradox else -12
            nom_fich='sky'
            lienzo.create_image(pixx(nodoj.x)+offsetx, ypx0, 
                        image=figs[nom_fich], anchor='n', tags=('apoyo'))
        if tramo.setj.bv_restri[4].get():
            centradox=True if (nodoj.bv_gdl[2].get() or itramo==ntramos-1) else False
            offsetx=0 if centradox else -12
            nom_fich='skt'
            lienzo.create_image(pixx(nodoj.x)+offsetx, ypx0, 
                        image=figs[nom_fich], anchor='n', tags=('apoyo'))



def do_nonulos():
    # especificar asientos & resortes:
    global entries_vistas
    
    if entries_vistas: # recoger datos, quitar entries y re-habilitar botones

        for itramo in range(ntramos):
            tramo=tramos[itramo]
            nodoi, nodoj = nodos[itramo], nodos[itramo+1]
            # recoger los del nodo i:
            for i in range(5):
                if tramo.seti.bv_restri[i].get():
                    try:
                        tramo.seti.valores[i]= float(tramo.seti.entries[i].get())
                    except (TypeError,ValueError):
                        texto='''No se entiende un dato de valor no nulo en un extremo izquierdo.'''
                        messagebox.showinfo(parent=v0, title='MAL', message=texto)
                        return
            # recoger los del nodo j:
            for i in range(5):
                if tramo.setj.bv_restri[i].get():
                    try:
                        tramo.setj.valores[i]= float(tramo.setj.entries[i].get())
                    except (TypeError,ValueError):
                        texto='''No se entiende un dato de valor no nulo en un extremo derecho.'''
                        messagebox.showinfo(parent=v0, title='MAL', message=texto)
                        return

        # comprobar que no se han puesto dos valores !=0 si gdl compartido
        ok=True
        for inodo in range(1,nnodos-1):
            nodo=nodos[inodo]
            tramo0, tramo1 = tramos[inodo-1], tramos[inodo]
            if nodo.bv_gdl[0].get():
                if tramo0.setj.bv_restri[0].get(): # and tramo1.seti.restri[0]
                    if tramo0.setj.valores[0] and tramo1.seti.valores[0]:
                        ok=False
                        break
            if nodo.bv_gdl[1].get():
                if tramo0.setj.bv_restri[1].get():
                    if tramo0.setj.valores[1] and tramo1.seti.valores[1]:
                        ok=False
                        break
                if tramo0.setj.bv_restri[3].get:
                    if tramo0.setj.valores[3] and tramo1.seti.valores[3]:
                        ok=False
                        break
            if nodo.bv_gdl[2].get():
                if tramo0.setj.bv_restri[2].get():
                    if tramo0.setj.valores[2] and tramo1.seti.valores[2]:
                        ok=False
                        break
                if tramo0.setj.bv_restri[4].get():
                    if tramo0.setj.valores[4] and tramo1.seti.valores[4]:
                        ok=False
                        break

        if not ok:
            texto='Hay dos valores no nulos asociados a un gdl compartido.'
            texto2='''Ponga en una de las casillas el valor no nulo del asiento o resorte (se tomará como el único representativo) y deje la otra con 0.0.'''
            messagebox.showinfo(parent=v0, title='MAL', message=texto,
                    detail=texto2)
            return


        frame_lienzo.delete('entry_s')
        btn_do_apoyos['state'] = 'normal'
        for tramo in tramos:
            tramo.seti.btn_ux.state(['!disabled'])
            tramo.seti.btn_uy.state(['!disabled'])
            tramo.seti.btn_fi.state(['!disabled'])
            tramo.seti.btn_ky.state(['!disabled'])
            tramo.seti.btn_kt.state(['!disabled'])
        
            tramo.setj.btn_ux.state(['!disabled'])
            tramo.setj.btn_uy.state(['!disabled'])
            tramo.setj.btn_fi.state(['!disabled'])
            tramo.setj.btn_ky.state(['!disabled'])
            tramo.setj.btn_kt.state(['!disabled'])
        

    else: # deshabilitar botones, poner entries y esperar nueva pulsacion
        
        btn_do_apoyos['state'] = 'disabled'
        for tramo in tramos:
            tramo.seti.btn_ux.state(['disabled'])
            tramo.seti.btn_uy.state(['disabled'])
            tramo.seti.btn_fi.state(['disabled'])
            tramo.seti.btn_ky.state(['disabled'])
            tramo.seti.btn_kt.state(['disabled'])
        
            tramo.setj.btn_ux.state(['disabled'])
            tramo.setj.btn_uy.state(['disabled'])
            tramo.setj.btn_fi.state(['disabled'])
            tramo.setj.btn_ky.state(['disabled'])
            tramo.setj.btn_kt.state(['disabled'])
        
        for itramo in range(ntramos):
            tramo=tramos[itramo]
            nodoi, nodoj = nodos[itramo], nodos[itramo+1]
            # el izdo (entries a la dcha):
            for i in range(5):
                if tramo.seti.bv_restri[i].get():
                    tramo.seti.grid_un_entry_s(i, nodoi.x, izda=False)
                    if tramo.seti.valores[i]:
                        tramo.seti.entries[i].delete(0,'end')
                        tramo.seti.entries[i].insert(0,str(tramo.seti.valores[i]))
            # el dcho (entries a la izda):
            for i in range(5):
                if tramo.setj.bv_restri[i].get():
                    tramo.setj.grid_un_entry_s(i, nodoj.x, izda=True)
                    if tramo.setj.valores[i]:
                        tramo.setj.entries[i].delete(0,'end')
                        tramo.setj.entries[i].insert(0,str(tramo.setj.valores[i]))

    # actualizar para proxima pulsacion
    entries_vistas = not entries_vistas


def pinta_cargasM(lienzo, ypx0, acota):

    max_fy=0.
    for tramo in tramos: # la fy
        a= max(abs(tramo.carga[3][0]),abs(tramo.carga[3][1]))
        if a>max_fy: max_fy=a
    if max_fy != 0:
        scal_fy=105/max_fy
        for tramo in tramos:
            if any(tramo.carga[3]):
                tramo.pinta_fy(lienzo, ypx0, scal_fy, acota)

    for tramo in tramos: # la mz
        if any(tramo.carga[4]):
            tramo.pinta_mz(lienzo,ypx0, acota)

    max_py= 0.
    for tramo in tramos: # la py
        a= max(abs(tramo.carga[1][0]),abs(tramo.carga[1][1]))
        if a> max_py: max_py=a 
    if max_py !=0:
        scal_py=130/max_py
        for tramo in tramos:
            if any((tramo.carga[1][0], tramo.carga[1][1])):
                tramo.pinta_py(lienzo, ypx0, scal_py, acota)


def pinta_cargasN(lienzo, ypx0, yoffset, acota):

    max_px=0.
    for tramo in tramos: # la px
        a= max(abs(tramo.carga[0][0]),abs(tramo.carga[0][1]))
        if a> max_px: max_px=a 
    if max_px !=0:
        scal_px= 80/max_px
        for tramo in tramos:
            if any((tramo.carga[0][0], tramo.carga[0][1])):
                tramo.pinta_px(lienzo, ypx0, scal_px, yoffset, acota)

    max_fx=0.
    for tramo in tramos: # la fx
        a= max(abs(tramo.carga[2][0]),abs(tramo.carga[2][1]))
        if a>max_fx: max_fx=a
    if max_fx != 0:
        scal_fx= 90/max_fx
        for tramo in tramos:
            if any(tramo.carga[2]):
                tramo.pinta_fx(lienzo, ypx0, scal_fx, acota)


def do_cargas():
    global tramos, nodos
    
    def limpiar_cargas():
        #nonlocal filas_carga
        for ifila in range(ntramos):
            for j in range(5):
                for k in range(2):
                    colus_carga[ifila][j][k].delete(0,'end')
                    colus_carga[ifila][j][k].insert(0, '0.0')
        hecho_cargas() # que guarde en los tramos y redibuje

    
    def hecho_cargas():
        global tramos, nodos

        # guarda las cargas en los tramos():

        for itramo in range(ntramos):
            tramo=tramos[itramo]
            tramo.carga=[[0.,0.],[0.,0.],[0.,0.],[0.,0.],[0.,0.]] # los tipos mutables
            for i in range(5): # px py fx fy mz
                for j in range(2): # izdo, dcho
                    try:
                        tramo.carga[i][j]=float(strv_carga[itramo][i][j].get())
                    except (TypeError, ValueError):
                        texto=f'No se entiende dato en el tramo {itramo}.'
                        messagebox.showinfo(parent=v2, title='MAL', message=texto)
                        return()

        # dibujar las cargas:
        
        frame_lienzo.delete('carga')
        pinta_cargasM(frame_lienzo, 300, acota=False)
        pinta_cargasN(frame_lienzo, 300, yoffset=100, acota=False)

        v2.destroy()
    
    
    # ventana de especificar cargas:
    
    v2=Toplevel(v0)
    v2.title('Cargas')

    frame_cargas= ttk.Frame(v2) # todas las cargas van a ser de tramo
    frame_cargas.grid(row=0, column=0, columnspan=2, padx=6, pady=6)
    texto ='El signo * indica que los gdl afectados están compartidos,\n'
    texto+='por lo que es indiferente a qué extremo de tramo se asigne\n'
    texto+='la carga. Se recomienda que, en su caso, ponga el valor de\n'
    texto+='dicha carga en cualquiera de las dos casillas adyacentes\n'
    texto+='dejando 0.0 en la otra. Si pone valores en ambas los mismos\n'
    texto+='se sumarán en el cálculo, aunque en la representación solo\n'
    texto+='se apreciará la mayor porque se dibujan centradas en el nodo.'
    CreateToolTip(frame_cargas, texto)
    
    
    ttk.Label(frame_cargas, text='tramo').grid(row=0, column=0)
    for i in range(ntramos):
        ttk.Label(frame_cargas, text=str(i)).grid(row=0, column=3*i+1, columnspan=3)
        ttk.Label(frame_cargas, text='iz').grid(row=1, column=3*i+1)
        ttk.Label(frame_cargas, text='dc').grid(row=1, column=3*i+2)
    
    ttk.Label(frame_cargas, text='px').grid(row=2, column=0)
    ttk.Label(frame_cargas, text='py').grid(row=3, column=0)
    ttk.Label(frame_cargas, text='Fx').grid(row=4, column=0)
    ttk.Label(frame_cargas, text='Fy').grid(row=5, column=0)
    ttk.Label(frame_cargas, text='Mz').grid(row=6, column=0)

    strv_carga=[]
    colus_carga=[]
    asteriscos=[]
    for itramo in range(ntramos):
        nodoj=nodos[itramo+1]
        colu_carga=[]       # corresponde a un tramo
        strv_carga_tramo=[]
        aster_tramo=[]
        for i in range(5):  # px py fx fy mz
            strv_carga_tramo.append([StringVar(value='0.0'), StringVar(value='0.0')])
            e1 = ttk.Entry(frame_cargas, width=5, 
                            textvariable=strv_carga_tramo[i][0])
            e2 = ttk.Entry(frame_cargas, width=5, 
                            textvariable=strv_carga_tramo[i][1])
            colu_carga.append([e1,e2])
            
            e1.grid(row=i+2, column=3*itramo+1)
            e2.grid(row=i+2, column=3*itramo+2)
            if itramo != ntramos-1:
                texto='*' if i>1 and nodoj.bv_gdl[i-2].get() else ''
                e = ttk.Label(frame_cargas, width=1, text=texto)
                aster_tramo.append(e)
                e.grid(row=i+2, column=3*itramo+3, padx=1)

        colus_carga.append(colu_carga)
        strv_carga.append(strv_carga_tramo)
        asteriscos.append(aster_tramo)

    # con lo anterior colus_carga & strv_carga quedan con la misma estructura
    # xxx.[itramo][itipo 1-5][dcho/izdo 1-2]
    # en ateriscos[itramo][3] se guardan labels que con '*' indican que el gdl esta
    # compartido en el nodo, y debe ponerse el izdo o el dcho =0. Para fx fy mz solo.
    # En ppio es solo una indicacion, no afecta al calculo.
    
    # Restituimos los valores previos, si los hay:
    for itramo in range(ntramos):
        tramo=tramos[itramo]
        for i in range(5):
            for j in range(2):
                a= str(tramo.carga[i][j]) or '0.0'
                strv_carga[itramo][i][j].set(a) 

    ttk.Button(v2, text='Hecho', command= hecho_cargas).grid(row=1, 
                column=1,padx=4, pady=2, sticky='w')
    ttk.Button(v2, text='Limpiar', command= limpiar_cargas).grid(row=1,
                column=0,padx=4, pady=2, sticky='e')



    v2.mainloop()




def pixx(x):
    return(round(50+x*scal_px))



def di_tramos():
    global tramos


    def hecho_tramos():
        # define tramos[] y nodos[]
        global tramos, nodos, ntramos, nnodos, scal_px, ltotal
        global btn_do_cargas, btn_do_apoyos, btn_do_nonulos
        
        # encontrar ntramos & asignar nnodos, xnodos, ltotal, scal_px:
        ntramos, ltotal, xnodos = 0, 0., [0.]
        for jtramo in range(hard_ntramos):
            try:
                ltotal += float(filas_tramos[0][jtramo].get())
                xnodos.append(ltotal)
                ntramos +=1
            except (TypeError, ValueError):
                break
        nnodos= ntramos +1  # sera=len(xnodos)
        scal_px=900./ltotal

        # asignar la lista de objetos tramo. Recordar, en filas_tramos
        # estan como L,E,I,A para cada tramo POR COLUMNAS.
        tramos= [] 
        for j in range(ntramos):
            tramo=Tramo(ctes=[0.,0.,0.,0.], inodos=[j, j+1]) 
                            # los putos tipos mutables quieren init explic
            for i in range(4):
                try:
                    tramo.ctes[i]= float(filas_tramos[i][j].get())
                except (TypeError, ValueError):
                    texto = f'Error en los datos:\nTramo {j}'
                    messagebox.showinfo(parent=v1,title='MAL', message=texto)
                    return()
            tramos.append(tramo)
        
        
        # creamos la lista de objetos nodo(). Creo el 1º & el ultimo aunque los
        # dejare deshabilitados en el gui.
        nodos=[]
        for j in range(nnodos):
            nodo=Nodo(j,xnodos[j])
            nodos.append(nodo)
        
        
        ###########################
        # pintamos el setup basico 
        ###########################
        
        frame_lienzo.delete("all")
        frame_lienzo.create_line(50,300,950,300, width=3) # la viga
        
        milb=ttk.Label(frame_general, text='\nux\nuy\nfi', font='Roman 13',background='white')
        frame_lienzo.create_window(20,15, anchor='n',window=milb) 
        
        for nodo in nodos:
            px=pixx(nodo.x)
            frame_lienzo.create_line(px,20,px,530,width=1,dash='2 4 6 4') # linea vert
            frame_lienzo.create_window(px, 20, anchor='n', window=nodo.cuadro)
        
        # deshabilitamos el 1º y el ultimo
        nodos[0].btn_shareux.state(['disabled'])
        nodos[0].btn_shareuy.state(['disabled'])
        nodos[0].btn_sharefi.state(['disabled'])
        nodos[-1].btn_shareux.state(['disabled'])
        nodos[-1].btn_shareuy.state(['disabled'])
        nodos[-1].btn_sharefi.state(['disabled'])
        
        
        milb=ttk.Label(frame_general, text='ux\nuy\nfi\nky\nkt', font='Roman 13',background='white')
        frame_lienzo.create_window(20,397, anchor='n',window=milb) 
        
        for i in range(ntramos): # ponemos los sets de s:
            tramo=tramos[i]
            inodo, jnodo = i, i+1

            px=pixx(xnodos[inodo])+10
            frame_lienzo.create_window(px, 400, anchor='n', window=tramo.seti.cuadro)
            for j in range(5):
                tramo.seti.bv_restri[j].set(False)
            tramo.seti.btn_ux.state(['disabled'])
            tramo.seti.btn_uy.state(['disabled'])
            tramo.seti.btn_fi.state(['disabled'])
            tramo.seti.btn_ky.state(['disabled'])
            tramo.seti.btn_kt.state(['disabled'])
            
            px=pixx(xnodos[jnodo])-10
            frame_lienzo.create_window(px, 400, anchor='n', window=tramo.setj.cuadro)
            for j in range(5):
                tramo.setj.bv_restri[j].set(False)
            tramo.setj.btn_ux.state(['disabled'])
            tramo.setj.btn_uy.state(['disabled'])
            tramo.setj.btn_fi.state(['disabled'])
            tramo.setj.btn_ky.state(['disabled'])
            tramo.setj.btn_kt.state(['disabled'])
            
            xmedio= (nodos[inodo].x + nodos[jnodo].x)/2.
            
            unlbl=ttk.Label(frame_general, text=str(i), style='jc_white.TLabel')
            frame_lienzo.create_window(pixx(xmedio), 90, anchor='n', window=unlbl)


        # botones tk, NO ttk, de "incorporar" (uniones y apoyos), 
        # de "valores no nulos" (solo en apoyos),
        # y de definir cargas
        
        btn_do_uniones= Button(frame_general, text='incorporar\nuniones',background='#D7ECFF',width=8, height=2, command=do_uniones)
        frame_lienzo.create_window(980, 35, anchor='nw', window=btn_do_uniones)

        btn_do_cargas= Button(frame_general, text='definir\ncargas',background='#FFD4D0',width=8, height=2, command=do_cargas, state=['disabled'])
        frame_lienzo.create_window(980, 190, anchor='nw', window=btn_do_cargas)

        btn_do_apoyos= Button(frame_general, text='incorporar\napoyos',background='#DAFFD7',width=8, height=2, command=do_apoyos, state=['disabled'])
        frame_lienzo.create_window(980, 400, anchor='nw', window=btn_do_apoyos)

        btn_do_nonulos= Button(frame_general, text='valores\nno nulos',background='#DAFFD7',width=8, height=2, command=do_nonulos, state=['disabled'])
        frame_lienzo.create_window(980, 460, anchor='nw', window=btn_do_nonulos)

        texto ='Pulse una vez para especificar valores de resortes\n'
        texto+='y de asientos de la sustentación. Pulse otra vez\n'
        texto+='para que se guarden dichos valores.\n'
        texto+='Cuando el gdl afectado esté compartido introduzca\n'
        texto+='el valor no nulo en una sola de las casillas.'
        CreateToolTip(btn_do_nonulos, texto)


        v1.destroy()


    v1=Toplevel(v0)
    v1.title('Tramos')
    
    texto= 'Introduzca los datos de los tramos comenzando por el izquierdo.\n'
    texto+='Deje en blanco los espacios sobrantes'
    ttk.Label(v1,text=texto).grid(row=0,column=0,columnspan=hard_ntramos, padx=3, pady=3)

    texto=''
    for i in range(hard_ntramos):
        texto = f'tramo {i}  '
        ttk.Label(v1,text=texto).grid(row=1, column=i+1, padx=3, pady=3)

    ttk.Label (v1,text='L:').grid(row=2, column=0, padx=3, pady=3)
    ttk.Label (v1,text='E:').grid(row=3, column=0, padx=3, pady=3)
    ttk.Label (v1,text='I:').grid(row=4, column=0, padx=3, pady=3)
    ttk.Label (v1,text='A:').grid(row=5, column=0, padx=3, pady=3)

    # organizo [cte, tramo] para que se parezca a la disposicion fisica
    # y relleno con valores previos si los hay

    filas_tramos=[]
    for i in range(4):
        fila=[]
        for j in range(hard_ntramos):
            e= ttk.Entry(v1, width=9)
            e.grid(row=i+2, column=j+1, padx=3, pady=3)
            fila.append(e)
        filas_tramos.append(fila)
    
    for i in range(4):
        for j in range(ntramos):
            a=tramos[j].ctes[i]
            filas_tramos[i][j].insert(0,str(a))
    
    ttk.Button(v1, text='Hecho', command=hecho_tramos).grid(column=1, row=6, columnspan=hard_ntramos, padx=3, pady=5)

    
    v1.mainloop()


def cerrar_graf():
    global v3,v4
    try:
        v3.destroy()
        v4.destroy()
    except:
        pass


def a_cargar():
    
    
    
    nfcompleto= filedialog.askopenfilename(parent=v0, title='Abrir archivo')



def a_guardar():
    
    
    nfcompleto= filedialog.asksaveasfilename(parent=v0, title='Nombre del archivo a guardar')
    



def imprimir():
    global frame_lienzoM, frame_lienzoN

    nomdir = filedialog.askdirectory(parent=v0, title= 'Carpeta de destino')
    
    salida_texto(carpeta=nomdir)
    
    patz1, patz2 = path.join(nomdir,'diagramas_flexion.ps'), path.join(nomdir,'diagramas_axil.ps')
    frame_lienzoM.postscript(file=patz1, height='1600', pagewidth='16c')
    frame_lienzoN.postscript(file=patz2, height='1100', pageheight='16c')
    
    texto1= 'Archivos guardados correctamente:'
    texto2= '- "diagramas_flexion.ps" \n- "diagramas_axil.ps" \n- "jcontinua.out" '
    texto2+='\nCarpeta: ' + nomdir
    messagebox.showinfo(message=texto1, detail=texto2, title='BIEN',parent=v0)




def salida_grafica():
    global tramos, v3,v4, frame_lienzoN, frame_lienzoM

    # En la salida grafica se vuelven a calcular las reacciones, esta vez por equilibrio 
    # directo. Es innecesario porque ya las dejo calcular() en tr.reacc[], pero lo hice para
    # comprobar que salía igual... y asi quedo.
    
    def reacM():
        
        def pintar_una_reacM(tipo, x, valor, centrado ):
            # tipo= V o M; x la del nodo (en px), valor=float, centrado= True o 'dcha' o 'izda'
            
            if centrado==True:
                xoffset=0
            elif centrado== 'dcha' :
                xoffset= 12
            else : # centrado es = izda
                xoffset=-12
            x0 = x+ xoffset

            if tipo=='V':
                yini, yfin = ypx0_M[0]+49, ypx0_M[0]+109
                if valor<0 : yini, yfin= yfin, yini
                frame_lienzoM.create_line(x0, yini, x0, yfin, arrow='last', 
                                            width=4, fill='#E342FF', arrowshape=(10,16,6))
                # acota:
                frame_lienzoM.create_text(x+4+xoffset, (yini+yfin)/2+rr(-9,9),
                    font=('Times','12','bold'), fill='#E342FF', 
                    text='{:.3g}'.format(abs(valor)), anchor='w')
            
            elif tipo=='M':
                signo= np.sign(valor)
                a, b = [x-25+xoffset, ypx0_M[0]+36], [x+25+xoffset, ypx0_M[0]+66]
                frame_lienzoM.create_arc(a[0],a[1],b[0],b[1], outline='#E342FF',
                            style='arc', start=184, extent=172, width=3)

                x0= b[0] if signo >0 else a[0]
                x1= b[0]+8 if signo >0 else a[0]-8
                y0, y1 = ypx0_M[0]+56, ypx0_M[0]+50
                frame_lienzoM.create_line(x0, y0, x1, y1, 
                        arrow='last', width=3, fill='#E342FF', arrowshape=(10,16,6))
                # acota:
                frame_lienzoM.create_text(b[0]+6, ypx0_M[0]+50+rr(-9,9),font=('Times','12','bold'), 
                        fill='#E342FF', anchor='w', text='{:.3g}'.format(abs(valor)))

        for inodo in range (1,nnodos-1):
            nodo = nodos[inodo]
            trs= [ tramos[inodo-1], tramos[inodo] ]

            # la fy:
            if nodo.bv_gdl[1].get(): # uy compartido
                centrado=True
                reac= (-trs[0].carga[3][1]- trs[1].carga[3][0]
                       +trs[0].fa_tr[1][0]+ trs[1].fa_tr[0][0] )
                       
                if trs[0].setj.bv_restri[1].get() or trs[1].seti.bv_restri[1].get(): # apoyo
                    pintar_una_reacM( 'V', pixx(nodo.x), reac, centrado)  

                elif trs[0].setj.bv_restri[3].get() or trs[1].seti.bv_restri[3].get(): # resorte
                    k= trs[0].setj.valores[3] or trs[1].seti.valores[3]
                    pintar_una_reacM( 'V', pixx(nodo.x), reac, centrado) 
                else:
                    pass
                    #print(f'\nNodo {inodo}. Equilibrio y (~0?) = {reac}') # no restringido

            else: # uy no compartido
                if trs[0].setj.bv_restri[1].get() or trs[0].setj.bv_restri[3].get():
                    centrado='izda'
                    reac= -trs[0].carga[3][1]+ trs[0].fa_tr[1][0]
                    pintar_una_reacM( 'V', pixx(nodo.x), reac, centrado) 
                if trs[1].seti.bv_restri[1].get() or trs[1].seti.bv_restri[3].get():
                    centrado='dcha'
                    reac= -trs[1].carga[3][0]+ trs[1].fa_tr[0][0]
                    pintar_una_reacM( 'V', pixx(nodo.x), reac, centrado) 

            # la mz:
            if nodo.bv_gdl[2].get(): # fi compartido
                centrado=True
                reac= (-trs[0].carga[4][1]- trs[1].carga[4][0]
                       +trs[0].fa_tr[1][1]+ trs[1].fa_tr[0][1] )

                if trs[0].setj.bv_restri[2].get() or trs[1].seti.bv_restri[2].get(): # apoyo
                    pintar_una_reacM( 'M', pixx(nodo.x), reac, centrado) 

                elif trs[0].setj.bv_restri[4].get() or trs[1].seti.bv_restri[4].get(): # resorte
                    k= trs[0].setj.valores[4] or trs[1].seti.valores[4]
                    pintar_una_reacM( 'M', pixx(nodo.x), reac, centrado) 
                else:
                    pass
                    #print(f'\nNodo {inodo}. Equilibrio M (~0?) = {reac}')

            else: # fi no compartido
                if trs[0].setj.bv_restri[2].get() or trs[0].setj.bv_restri[4].get():
                    centrado='izda'
                    reac= -trs[0].carga[4][1]+ trs[0].fa_tr[1][1]
                    pintar_una_reacM( 'M', pixx(nodo.x), reac, centrado) 
                if trs[1].seti.bv_restri[2].get() or trs[1].seti.bv_restri[4].get():
                    centrado='dcha'
                    reac= -trs[1].carga[4][0]+ trs[1].fa_tr[0][1]
                    pintar_una_reacM( 'M', pixx(nodo.x), reac, centrado) 
        
        
        # nodo 0, fy:
        tramo=tramos[0]
        if tramo.seti.bv_restri[1].get() or tramo.seti.bv_restri[3].get():
            reac= -tramo.carga[3][0] + tramo.fa_tr[0][0] 
            pintar_una_reacM( 'V', pixx(nodos[0].x), reac, True) 
        # nodo0, mz
        if tramo.seti.bv_restri[2].get() or tramo.seti.bv_restri[4].get():
            reac= -tramo.carga[4][0] + tramo.fa_tr[0][1] 
            pintar_una_reacM( 'M', pixx(nodos[0].x), reac, True) 

        # nodo n, fy:
        tramo=tramos[ntramos-1]
        if tramo.setj.bv_restri[1].get() or tramo.setj.bv_restri[3].get():
            reac= -tramo.carga[3][1] + tramo.fa_tr[1][0] 
            pintar_una_reacM( 'V', pixx(nodos[nnodos-1].x), reac, True)
        # nodo n, mz:
        if tramo.setj.bv_restri[2].get() or tramo.setj.bv_restri[4].get():
            reac= -tramo.carga[4][1] + tramo.fa_tr[1][1] 
            pintar_una_reacM( 'M', pixx(nodos[nnodos-1].x), reac, True)
        


    def reacN():
        
        def pintar_una_reacN( xn, valor, centrado):
            # x del nodo en px. reac float. centrado str o bool
            
            if centrado==True:
                xoffset=0
            elif centrado== 'dcha' :
                xoffset= 12
            else : # centrado es = izda
                xoffset=-12

            xini = xn+ xoffset
            xfin = xini+44 if reac>0 else xini-44
            if np.isclose(xn,50.) and valor<0 : xini, xfin= xn+44, xn

            if centrado=='dcha': # decalajes verticales para legibilidad
                y0 = ypx0_N[0] +90+20
            elif centrado=='izda':
                y0 = ypx0_N[0] +90-20
            else:
                y0 = ypx0_N[0] +90
            
            frame_lienzoN.create_line(xini,y0, xfin,y0, arrow='last', 
                                        width=4, fill='#E342FF', arrowshape=(10,16,6))
            # acota:
            frame_lienzoN.create_text((xini+xfin)/2, y0-14, font=('Times','12','bold'),
                            fill='#E342FF', text='{:.3g}'.format(abs(valor)))


        for inodo in range (1, nnodos-1):
            nodo = nodos[inodo]
            trs= [ tramos[inodo-1], tramos[inodo] ]

            # dibujo las fx de los nodos intermedios desde el apoyo:

            if nodo.bv_gdl[0].get(): # ux compartido
                centrado=True
                reac= (-trs[0].carga[2][1]- trs[1].carga[2][0]
                       +trs[0].fa_tr[1][4]+ trs[1].fa_tr[0][4] )
                if trs[0].setj.bv_restri[0].get() or trs[1].seti.bv_restri[0].get(): # apoyo
                    pintar_una_reacN( pixx(nodo.x), reac, centrado)  
                            # (no hay resortesN)
                else:
                    pass
                    #print(f'\nNodo {inodo}. Equilibrio x (~0?) = {reac}') # no restringido

            else: # ux no compartido
                if trs[0].setj.bv_restri[0].get() :
                    centrado='izda'
                    reac= -trs[0].carga[2][1]+ trs[0].fa_tr[1][4]
                    pintar_una_reacN( pixx(nodo.x), reac, centrado) 
                if trs[1].seti.bv_restri[0].get() :
                    centrado='dcha'
                    reac= -trs[1].carga[2][0]+ trs[1].fa_tr[0][4]
                    pintar_una_reacN( pixx(nodo.x), reac, centrado) 
        
        # nodo 0, fx:
        tramo=tramos[0]
        if tramo.seti.bv_restri[0].get() :
            reac= -tramo.carga[2][0] + tramo.fa_tr[0][4] 
            pintar_una_reacN( pixx(nodos[0].x), reac, True) 
        
        # nodo n, fx:
        tramo=tramos[ntramos-1]
        if tramo.setj.bv_restri[0].get() :
            reac= -tramo.carga[2][1] + tramo.fa_tr[1][4] 
            pintar_una_reacN( pixx(nodos[nnodos-1].x), reac, True)
    
    
    


    def trazar_una_grafica (igr):  # igr 0->4 : V, M, fi, u ;  ventana v3
        
        lienzo = frame_lienzoM if igr<4 else frame_lienzoN
        ypx0= ypx0_M[igr+1]    if igr<4 else ypx0_N[igr-3]
        
        lienzo.create_line(50,ypx0,950,ypx0, width=2)

        if maximos[igr]==0. :
            #print('\n Hallada grafica identicamente nula (no se dibuja).\n')
            return(1)
            
        scal_py= 145/maximos[igr]
        pixy= lambda y: round(ypx0+y*scal_py) 
        pixx= lambda x: 50 + (x0+x)*scal_px # redefino pixx (scal_px es el global)

        for itramo in range(ntramos):
            tramo=tramos[itramo]
            funcion=tramo.funciones[igr]
            x0= nodos[tramo.inodos[0]].x
            
            # trazar el dibujo basico:
            lista_pa_trazar=[pixx(0),pixy(0.)]
            x_val, y_val= [],[]
            for x in np.linspace(0.,tramo.ctes[0],20): # 12 o lo que vaya bien
                y=funcion(x)
                x_val.append(x)
                y_val.append(y)
                lista_pa_trazar.append(pixx(x))
                lista_pa_trazar.append(pixy(y))
            lista_pa_trazar.append(pixx(tramo.ctes[0]))
            lista_pa_trazar.append(pixy(0.))

            lienzo.create_polygon(lista_pa_trazar, fill='#FFFEC3', outline='black')

            # acotar extremos de tr. Cifras siempre dentro del tr (que lo amarillo no tape)
            
            # acotar extremo izdo del tr  (salvo =0 ó que el anterior coincida)
            poco= maximos[igr]/1.e5
            if itramo == 0: 
                a=tramo.fa_tr[0][igr]
                if abs(a) > poco:
                    c='sw' if lista_pa_trazar[3]<0. else'nw'
                    lienzo.create_text(lista_pa_trazar[2], lista_pa_trazar[3],
                                        text='{:.3g}'.format(abs(a)), anchor=c)
            else: # comprobar que no coincide con el anterior (con -el si es VMN)
                a,b = tramo.fa_tr[0][igr], tramos[itramo-1].fa_tr[1][igr]
                if igr in (0,1,4) : b=-b
                if not np.isclose(a,b) and (abs(a) > poco):
                    c='sw' if lista_pa_trazar[3]<0. else'nw'
                    lienzo.create_text(lista_pa_trazar[2], lista_pa_trazar[3],
                                        text='{:.3g}'.format(abs(a)), anchor=c)
            # acotar extremo dcho del tr  (salvo =0; ya coincida el posterior o no)
            a=tramo.fa_tr[1][igr] 
            if abs(a) > poco:
                c='se' if lista_pa_trazar[-3]<0. else'ne'
                lienzo.create_text(lista_pa_trazar[-4], lista_pa_trazar[-3],
                                    text='{:.3g}'.format(abs(a)), anchor=c)

            # acotar maximos si los hay
            
            haymax, xdemax = False, []
            if igr==0: # el del cortante hay que calcularlo ahora:
                if ((tramo.carga[1][0]>0. and tramo.carga[1][1]<0) or
                    (tramo.carga[1][0]<0. and tramo.carga[1][1]>0)):
                    haymax=True
                    xdemax.append(  abs(tramo.carga[1][0])*tramo.ctes[0]/
                              (abs(tramo.carga[1][0])+abs(tramo.carga[1][1]) )  )
            elif igr==4: # el de N tambien
                if ((tramo.carga[0][0]>0. and tramo.carga[0][1]<0) or
                    (tramo.carga[0][0]<0. and tramo.carga[0][1]>0)):
                    haymax=True
                    xdemax.append(  abs(tramo.carga[0][0])*tramo.ctes[0]/
                              (abs(tramo.carga[0][0])+abs(tramo.carga[0][1]) )  )
            else:
                haymax=bool(tramo.haycero[igr-1])
                xdemax=tramo.haycero[igr-1] # puede haber varios, es una lista

            if haymax:
                for xm in xdemax:
                    if not(np.isclose(xm, tramo.ctes[0])) and xm>tramo.ctes[0]/1.e3:
                        ym= tramo.funciones[igr](xm)
                        c='s' if ym<0. else'n'
                        lienzo.create_text( pixx(xm), pixy(ym),
                                                text='{:.3g}'.format(abs(ym)), anchor=c)
                        c='n' if ym<0. else's'
                        lienzo.create_text( pixx(xm), pixy(0.),
                                                text='   x=\n{:.3g}'.format(x0+xm), anchor=c)
                        lienzo.create_line( pixx(xm), pixy(0.), pixx(xm), pixy(ym),
                                            width=1, fill='#C622CC', dash=(6,2,2,2) )

            # poner simbolitos de + -

            letra, L = ['v','m','g','', 'n','w'], tramo.ctes[0] 
            if igr !=3: # para el despl. no hay simbolo
                # poner simbolos en todos los trozos entre ceros me parece excesivo, y ademas
                # va a ser raro que en fi (o anterior) haya mas de un cero. Opto por:
                # - si no hay cero, poner un simbolo.
                # - si hay 1 cero, poner dos simbolos
                # - si hay mas de un cero, poner un simbolo entre el 1º y el 2º 

                if len(tramo.haycero[igr])==1:
                    xdecero= tramo.haycero[igr][0]
                    xa, xb = xdecero/2,  (L+xdecero)/2
                    ya, yb= tramo.funciones[igr](xa), tramo.funciones[igr](xb)
                    c='p' if ya>0. else 'n'
                    nomfig='d'+letra[igr]+c
                    lienzo.create_image(pixx(xa), pixy(ya/2), image= figs[nomfig])
                    c='p' if yb>0. else 'n'
                    nomfig='d'+letra[igr]+c
                    lienzo.create_image(pixx(xb), pixy(yb/2), image= figs[nomfig])
                elif len(tramo.haycero[igr])==0: # fig a 1/3 del extremo mayor
                    xa= L/3 if abs(tramo.fa_tr[0][igr])>abs(tramo.fa_tr[1][igr]) else 2*L/3
                    ya= tramo.funciones[igr](xa)
                    c='p' if ya>0. else 'n'
                    nomfig='d'+letra[igr]+c
                    lienzo.create_image(pixx(xa), pixy(ya/2), image= figs[nomfig])
                else: # hay mas de un cero
                    xa=(tramo.haycero[igr][0]+tramo.haycero[igr][1])/2
                    ya = tramo.funciones[igr](xa)
                    c='p' if ya>0. else 'n'
                    nomfig='d'+letra[igr]+c
                    lienzo.create_image(pixx(xa), pixy(ya/2), image= figs[nomfig])



    # preparar las ventanas:
    
    v3 = Toplevel(v0)
    v3.title('Flexión')
    frame_lienzoM= Canvas(v3, width=1100, height=700, background='white', scrollregion=(0,0,1000,1600))
    vbarM=Scrollbar(v3,orient=VERTICAL)
    vbarM.grid(column=1, row=0, sticky='ns')
    vbarM.config(command=frame_lienzoM.yview)
    frame_lienzoM.config(width=1100,height=700)
    frame_lienzoM.config(yscrollcommand=vbarM.set)
    frame_lienzoM.grid(column=0, row=0, sticky='nsew')


    v4 = Toplevel(v0)
    v4.title('Tracc-Compr')
    frame_lienzoN= Canvas(v4, width=1100, height=700, background='white', scrollregion=(0,0,1000,1100))
    vbarN=Scrollbar(v4,orient=VERTICAL)
    vbarN.grid(column=1, row=0, sticky='ns')
    vbarN.config(command=frame_lienzoN.yview)
    frame_lienzoN.config(width=1100,height=700)
    frame_lienzoN.config(yscrollcommand=vbarN.set)
    frame_lienzoN.grid(column=0, row=0, sticky='nsew')


    # encontrar extremos groseros para dibujar. Muestreo en 0, L/4, L/2, 3L/4, L
    
    extremos=[  [ 1.e99, 1.e99, 1.e99, 1.e99, 1.e99, 1.e99],
                [-1.e99,-1.e99,-1.e99,-1.e99,-1.e99,-1.e99]  ] # minimos, maximos
    maximos = [0., 0., 0., 0., 0., 0.] # en valor abs, para scal
    
    for tramo in tramos:
        for i in range(6): # para V M fi u N w
            b= [ a for a in tramo.funciones[i](np.linspace(0,tramo.ctes[0],5)) ]
            if max(b) > extremos[1][i]: extremos[1][i]=max(b)
            if min(b) < extremos[0][i]: extremos[0][i]=min(b)
    maximos= [ max(abs(extremos[0][i]), abs(extremos[1][i])) for i in range(6) ]


    # Asignar cotas y-pixel base para las graficas. La viga ocupara 0-150-300. Las graficas, si
    # el posit_ o el neg_ no tiene nada le asigno 20: -> 0-20-170 ó 0-150-170 
    # si hay ambos asigno 150 al mayor & proporcianal +10 al menor
    
    ypx0_M, ypx0_N =[150], [150] 
    ocupa_M, ocupa_N =[290],[290]

    for igr in range(4): # para las V M fi u
        if extremos[0][igr]==0 and extremos[1][igr]==0:
            ypx0_M.append(sum(ocupa_M)+20)
            ocupa_M.append(40)
        elif extremos[1][igr] <=0: # todo es zona -
            ypx0_M.append(sum(ocupa_M)+150)
            ocupa_M.append(180)
        elif  extremos[0][igr] >= 0: # todo es zona +
            ypx0_M.append(sum(ocupa_M)+20)
            ocupa_M.append(180)
        else: # hay ambos hemisferios

            if abs(extremos[0][igr])>abs(extremos[1][igr]): # el negativo es mayor
                ypx0_M.append(sum(ocupa_M)+150)
                c= 10+round(150*abs(extremos[1][igr])/abs(extremos[0][igr])) # pixs del posit 
                ocupa_M.append(170+c)
            else: # el positivo es mayor
                c= 10+round(170*abs(extremos[0][igr])/abs(extremos[1][igr])) # pixs del negat 
                ypx0_M.append(sum(ocupa_M)+c)
                ocupa_M.append(170+c)

    for igr in range(2): # para las N, w
        if extremos[0][igr+4]==0 and extremos[1][igr+4]==0:
            ypx0_N.append(sum(ocupa_N)+20)
            ocupa_N.append(40)
        elif extremos[1][igr+4] <=0: # todo es zona -
            ypx0_N.append(sum(ocupa_N)+150)
            ocupa_N.append(180)
        elif  extremos[0][igr+4] >= 0: # todo es zona +
            ypx0_N.append(sum(ocupa_N)+20)
            ocupa_N.append(180)
        else: # hay ambos hemisferios

            if abs(extremos[0][igr+4])>abs(extremos[1][igr+4]): # el negativo es mayor
                ypx0_N.append(sum(ocupa_N)+150)
                c= 12+round(150*abs(extremos[1][igr+4])/abs(extremos[0][igr+4])) 
                ocupa_N.append(170+c)
            else: # el positivo es mayor
                c= 12+round(150*abs(extremos[0][igr+4])/abs(extremos[1][igr+4])) # pixs del posit
                ypx0_N.append(sum(ocupa_N)+c)
                ocupa_N.append(170+c)


    # reproducir el dibujo de cargas & cc flex en v3
    frame_lienzoM.create_line(50,ypx0_M[0],950,ypx0_M[0], width=3) # la viga
    pinta_uniones ( frame_lienzoM, ypx0_M[0] )
    pinta_apoyos  ( frame_lienzoM, ypx0_M[0] )
    pinta_cargasM ( frame_lienzoM, ypx0_M[0], acota=True )

    # reproducir el dibujo de cargas & cc axil en v4
    frame_lienzoN.create_line(50,ypx0_N[0],950,ypx0_N[0], width=3) # la viga
    pinta_uniones ( frame_lienzoN, ypx0_N[0] )
    pinta_apoyos  ( frame_lienzoN, ypx0_N[0] )
    pinta_cargasN ( frame_lienzoN, ypx0_N[0], yoffset=0, acota=True )
    
    # trazar lineas verticales
    for nodo in nodos:
        x= pixx(nodo.x)
        frame_lienzoM.create_line(x, ypx0_M[0]+90, x, ypx0_M[4]+60, width=2, dash='2 4 6 4') 
        frame_lienzoN.create_line(x, ypx0_N[0]+90, x, ypx0_N[2]+60, width=2, dash='2 4 6 4') 
    
    for i in range(6): # graficas de V-M-fi-u -N-w
        trazar_una_grafica(i)
    
    reacM()
    reacN()

    return()


def salida_texto(carpeta=''):
    # Salida de texto, ya sea a terminal (por defecto) o a fichero jcontinua.out:
    if carpeta:
        nfcompleto= path.join(carpeta,'jcontinua.out')
        f=open(nfcompleto, 'w')
    else:
        f=sys.stdout

    nombrecillos=['V ','M ','fi','u ','N ','w ']
    print('#'*80,file=f)
    print('\nValores de las fuerzas, giros y desplazamientos en los extremos de tramo.',file=f)
    print(' V= fuerza vertical (positiva hacia abajo)',file=f)
    print(' M= momento en z    (positivo saliente)',file=f)
    print(' fi= giro en z      (positivo saliente)',file=f)
    print(' u= desplaz.vertic. (positivo hacia abajo)',file=f)
    print(' N= fuerza horiz.   (positiva hacia la derecha)',file=f)
    print(' u= desplaz.horiz.  (positivo hacia la derecha)',file=f)
    print('   Nota 1: en este listado se dan FUERZAS, no esfuerzos.',file=f)
    print('   Nota 2: en este listado las x_max son locales del tramo.',file=f)
    for itramo in range(ntramos):
        tramo=tramos[itramo]
        print(f'\nTramo: {itramo} ','-'*50,file=f)
        print(' '*14 + 'extr_i        extr_j    (x_max, valor)...',file=f)
        for igr in range(6):
            texto  = '    '+nombrecillos[igr]
            texto += '  {:12.5g}'.format(tramo.fa_tr[0][igr])
            texto += '  {:12.5g}'.format(tramo.fa_tr[1][igr])
            for c in tramo.maximos[igr]:
                texto += '    ({:.5g}, {:.5g})'.format(c[0],c[1])
            print(texto,file=f)

        print('  Reacciones (=0 si no hay ó si figura en tramo contiguo):',file=f) 
        print('    V ',' {:12.5g}  {:12.5g}'.format(tramo.reacc[0][0], tramo.reacc[1][0] ),
                        file=f)
        print('    M ',' {:12.5g}  {:12.5g}'.format(tramo.reacc[0][1], tramo.reacc[1][1] ),
                        file=f)
        print('    N ',' {:12.5g}  {:12.5g}'.format(tramo.reacc[0][2], tramo.reacc[1][2] ),
                        file=f)

    print()
    print('#'*80,file=f)
    
    if carpeta: f.close()



def calcular():
    # usa nodos[], tramos[], nnodos, ntramos, pero no los cambia (??)***
    
    def comprueba_f_apoyo():
        # que no haya fuerza concentrada en apoyo. Lo hago aqui porque antes
        # no se sabe si las restri[] & las carga[] son definitivas
        
        for itramo in range(ntramos):
            tramo=tramos[itramo]
            nodoi, nodoj = nodos[itramo],nodos[itramo+1]
            if (  ((tramo.seti.bv_restri[0].get() and tramo.carga[2][0]) or
                   (tramo.seti.bv_restri[1].get() and tramo.carga[3][0]) or
                   (tramo.seti.bv_restri[2].get() and tramo.carga[4][0])  )
                                        or
                  ((tramo.setj.bv_restri[0].get() and tramo.carga[2][1]) or
                   (tramo.setj.bv_restri[1].get() and tramo.carga[3][1]) or
                   (tramo.setj.bv_restri[2].get() and tramo.carga[4][1])  )   ):
                texto = 'Hay una carga concentrada sobre un apoyo. '
                texto2= 'Esto tiene poco sentido en el modelo de barras '
                texto2+='y el programa no lo admite. En casos reales puede ocurrir, '
                texto2+='pero sus efectos serán de tipo local (abolladura '
                texto2+='del alma, resistencia del apoyo, etc).'
                messagebox.showinfo(parent=v0, message=texto, detail=texto2, 
                        title='MAL')
                return(1)
        return(0)
    
    def mecanismo():
        # comprueba si es un mecanismo (cuenta global basica):
        gl= 3*ntramos
        for i in range(1,nnodos-1): # no el primero ni el ultimo
            nodo=nodos[i]
            for j in range(3):
                if nodo.bv_gdl[j].get(): 
                    gl -= 1

            tramo0, tramo1 = tramos[i-1], tramos[i]
            c=0
            for j in range(5):
                if tramo0.setj.bv_restri[j].get(): 
                    c+= 1
                if tramo1.seti.bv_restri[j].get(): 
                    if j==0 and not nodo.bv_gdl[0].get(): c += 1
                    if j==1 and not nodo.bv_gdl[1].get(): c += 1
                    if j==2 and not nodo.bv_gdl[2].get(): c += 1
                    if j==3 and not nodo.bv_gdl[1].get(): c += 1
                    if j==4 and not nodo.bv_gdl[2].get(): c += 1
            gl -= c
        
        for j in range(5): # para el primer y ultimo nodos
            if tramos[0].seti.bv_restri[j].get():         gl -= 1 
            if tramos[ntramos-1].setj.bv_restri[j].get(): gl -= 1
        
        mecanismo=True if gl>0 else False
        if mecanismo:
            texto = f'El sistema globalmente es un mecanismo (gdl={gl}).'
            texto2= 'Revise la sustentación y las uniones para asegurar '
            texto2+='que no lo sea. Cuide también que no haya mecanismos '
            texto2+='locales (el programa no lo comprueba explícitamente).'
            messagebox.showinfo(parent=v0, message=texto, detail=texto2, 
                        title='MAL')
        return(mecanismo)


    # En las funciones siguientes 0-1 indica el tramo anterior o el posterior
    # al nodo, mientras que i-j indica si en el tramo considerado es nodo i ó j.
    # Por eso solo hay 0j & 1i. La notación es redundante pero la mantengo porque
    # me ayuda a no equivocarme con los términos de mis tablas.
    # Las funciones en si mismas ensamblan en las matrices globales los esfuerzos
    # V ó M de extremos de barra.
    
    # Dado inodo tenemos:           itramo0=inodo-1 ;   itramo1=inodo
    # Sus ecs de equil estan en:    4*(inodo-1)+2   ; idem+3, +4, & +5
    # ordenadas como:
    #          si es no-compartido: fy0, mz0, fy1, mz1  
    #          si es compartido:    fy0, mz0,  u0j=u1i,  fi0j=fi1i
    
    # Las incog. de un tramo en:    4*itramo    ; idem+1, +2, & +3
    # ordenadas como:               ui, fii, uj, fij

    # iec = nº de ecuacion del K_glob donde sumaremos
    # en V0j & M0j siempre iec=0-1; V1i & M1i iran a ecs 0-1 si es compartido
    
    def pon_V0j(inodo, iec): # pone -V0j
        nonlocal K_glob, f_glob
        
        itramo=inodo-1
        tramo=tramos[itramo]
        iin= 4*itramo
        EI,L= tramo.ctes[1]*tramo.ctes[2], tramo.ctes[0]
        p,q = tramo.carga[1]

        a12, a6 = 12*EI/L**3, 6*EI/L**2
        K_glob[iec,iin+0] +=  a12
        K_glob[iec,iin+1] += -a6
        K_glob[iec,iin+2] += -a12
        K_glob[iec,iin+3] += -a6
        f_glob[iec] += -0.15*p*L -0.35*q*L
    
    def pon_M0j(inodo, iec): # pone -M0j
        nonlocal K_glob, f_glob
        
        itramo=inodo-1
        tramo=tramos[itramo]
        iin= 4*itramo
        EI,L= tramo.ctes[1]*tramo.ctes[2], tramo.ctes[0]
        p,q = tramo.carga[1]

        a2, a4, a6 = 2*EI/L, 4*EI/L, 6*EI/L**2
        K_glob[iec,iin+0] +=  a6
        K_glob[iec,iin+1] += -a2
        K_glob[iec,iin+2] += -a6
        K_glob[iec,iin+3] += -a4
        f_glob[iec] += -q*L**2/20 - p*L**2/30        
        
    def pon_V1i(inodo, iec): # pone -V1i
        nonlocal K_glob, f_glob
        
        itramo=inodo
        tramo=tramos[itramo]
        iin= 4*itramo
        EI,L= tramo.ctes[1]*tramo.ctes[2], tramo.ctes[0]
        p,q = tramo.carga[1]
        
        a12, a6 = 12*EI/L**3, 6*EI/L**2
        K_glob[iec,iin+0] += -a12
        K_glob[iec,iin+1] +=  a6
        K_glob[iec,iin+2] +=  a12
        K_glob[iec,iin+3] +=  a6
        f_glob[iec] += -0.15*q*L -0.35*p*L
        
    def pon_M1i(inodo, iec): # pone -M1i
        nonlocal K_glob, f_glob
        
        itramo=inodo
        tramo=tramos[itramo]
        iin= 4*itramo
        EI,L= tramo.ctes[1]*tramo.ctes[2], tramo.ctes[0]
        p,q = tramo.carga[1]

        a2, a4, a6 = 2*EI/L, 4*EI/L, 6*EI/L**2
        K_glob[iec,iin+0] +=  a6
        K_glob[iec,iin+1] += -a4
        K_glob[iec,iin+2] += -a6
        K_glob[iec,iin+3] += -a2
        f_glob[iec] +=  p*L**2/20 + q*L**2/30    

    cerrar_graf()

    # comprobaciones previas:
    if comprueba_f_apoyo() or len(nodos)==0 or len(tramos)==0: return(1)
    if mecanismo() : return(1)
    

    # Comienza el calculo de flexion. 
    # Hay 4 incog por tramo: (uy fi)_i (uy fi)_j
    K_glob= np.zeros((4*ntramos, 4*ntramos))
    a_glob= np.zeros(4*ntramos)
    f_glob= np.zeros(4*ntramos)
    prescindir= np.array([False]*(4*ntramos))

    # Procedo por nodos, equilibrando cada uno. El primero y el ultimo van aparte.
    for inodo in range (1,nnodos-1):
        nodo=nodos[inodo]
        iec= 4*(inodo-1)+2
        itramo0, itramo1= inodo-1, inodo
        iin0, iin1 = 4*itramo0, 4*itramo1
        
        if nodo.bv_gdl[1].get(): # gdl 1 compartido
            # ecuacion 1&3:
            pon_V0j(inodo, iec)
            pon_V1i(inodo, iec)
            ky= tramos[itramo0].setj.valores[3] or tramos[itramo1].seti.valores[3]
            Fext= tramos[itramo0].carga[3][1]+tramos[itramo1].carga[3][0]
            K_glob[iec,iin0+2] += -ky
            f_glob[iec] += -Fext

            K_glob[iec+2,iin0+2]=  1.
            K_glob[iec+2,iin1]  = -1.

        else: # gdl 1 no compartido
            # ecuacion 1:
            pon_V0j(inodo,iec)
            ky= tramos[itramo0].setj.valores[3]
            Fext= tramos[itramo0].carga[3][1]
            K_glob[iec,iin0+2] += -ky
            f_glob[iec] += -Fext
            # ecuacion 3:
            pon_V1i(inodo,iec+2)
            ky= tramos[itramo1].seti.valores[3]
            Fext= tramos[itramo1].carga[3][0]
            K_glob[iec+2,iin1] += -ky
            f_glob[iec+2] += -Fext

        if nodo.bv_gdl[2].get(): # gdl 2 compartido
            # ecuacion 2&4:
            pon_M0j(inodo, iec+1)
            pon_M1i(inodo, iec+1)
            kt= tramos[itramo0].setj.valores[4] or tramos[itramo1].seti.valores[4]
            Mext= tramos[itramo0].carga[4][1]+tramos[itramo1].carga[4][0]
            K_glob[iec+1, iin0+3] += -kt
            f_glob[iec+1] += -Mext

            K_glob[iec+3,iin0+3]=  1.
            K_glob[iec+3,iin1+1]= -1.

        else: # gdl 2 no compartido
            # ecuacion 2:
            pon_M0j(inodo,iec+1)
            kt= tramos[itramo0].setj.valores[4]
            Mext= tramos[itramo0].carga[4][1]
            K_glob[iec+1,iin0+3] += -kt
            f_glob[iec+1] += -Mext
            # ecuacion 4:
            pon_M1i(inodo,iec+3)
            kt= tramos[itramo1].seti.valores[4]
            Mext= tramos[itramo1].carga[4][0]
            K_glob[iec+3,iin1+1] += -kt
            f_glob[iec+3] += -Mext
            
    # El primer nodo:
    iec, iin = 0, 0
        # ecuacion 1:
    pon_V1i(0,iec)
    ky= tramos[0].seti.valores[3]
    Fext= tramos[0].carga[3][0]
    K_glob[iec,iin] += -ky
    f_glob[iec] += -Fext
        # ecuacion 2:
    pon_M1i(0,iec+1)
    kt= tramos[0].seti.valores[4]
    Mext= tramos[0].carga[4][0]
    K_glob[iec+1,iin+1] += -kt
    f_glob[iec+1] += -Mext
    

    # El ultimo nodo:
    iec, iin = 4*(nnodos-2)+2, 4*(ntramos-1)+2
        # ecuacion 1:
    pon_V0j(nnodos-1,iec)
    ky= tramos[ntramos-1].setj.valores[3]
    Fext= tramos[ntramos-1].carga[3][1]
    K_glob[iec,iin] += -ky
    f_glob[iec] += -Fext
        # ecuacion 2:
    pon_M0j(nnodos-1,iec+1)
    kt= tramos[ntramos-1].setj.valores[4]
    Mext= tramos[ntramos-1].carga[4][1]
    K_glob[iec+1,iin+1] += -kt
    f_glob[iec+1] +=  -Mext
        
    '''
    print('\nK & f globales:\n')
    for fila in K_glob:
        for k in fila:
            print('{:10.3g} '.format(k), end=' ')
        print('\n')

    print('\n')
    for f in f_glob:
        print ('{:10.3g} '.format(f), end=' ')
    print()
    '''
    
    # pongo ecuaciones prescindibles (& su gdl prescrito):

    for inodo in range(1, nnodos-1):
        nodo=nodos[inodo]
        tramo0, tramo1= tramos[inodo-1], tramos[inodo]
        iec= 4*(inodo-1) + 2

        if tramo0.setj.bv_restri[1].get(): # uy prescrito
            if nodo.bv_gdl[1].get(): # compartido
                valor= tramo0.setj.valores[1] or tramo1.seti.valores[1]
                a_glob[iec] = valor
                prescindir[iec] = True
            else: # no compartido
                valor= tramo0.setj.valores[1]
                a_glob[iec] = valor
                prescindir[iec] = True
                
        if tramo1.seti.bv_restri[1].get(): # uy prescrito
            if nodo.bv_gdl[1].get(): # esta compartido
                pass # dejar la ec de 1, -1 como esta (que calcule)
            else: # no esta compartido
                valor = tramo1.seti.valores[1]
                a_glob[iec+2] = valor
                prescindir[iec+2] = True
        
        if tramo0.setj.bv_restri[2].get(): # fi prescrito
            if nodo.bv_gdl[2].get(): # compartido:
                valor= tramo0.setj.valores[2] or tramo1.seti.valores[2]
                a_glob[iec+1] = valor
                prescindir[iec+1] = True
            else: # no compartido
                valor= tramo0.setj.valores[2]
                a_glob[iec+1] = valor
                prescindir[iec+1] = True
                
        if tramo1.seti.bv_restri[2].get(): # fi prescrito
            if nodo.bv_gdl[2].get(): # esta compartido
                pass # dejar la ec de 1, -1 como esta (que calcule)
            else: # no esta compartido
                valor = tramo1.seti.valores[2]
                a_glob[iec+3] = valor
                prescindir[iec+3] = True

    # el primer nodo:
    inodo, itramo, iec = 0, 0, 0
    nodo, tramo = nodos[inodo], tramos[itramo]
    if tramo.seti.bv_restri[1].get(): # uy prescrito
        valor= tramo.seti.valores[1]
        a_glob[iec] = valor
        prescindir[iec] = True
    if tramo.seti.bv_restri[2].get(): # fi prescrito
        valor= tramo.seti.valores[2]
        a_glob[iec+1] = valor
        prescindir[iec+1] = True

    # el ultimo nodo:
    inodo, itramo, iec = nnodos-1, ntramos-1, 4*(nnodos-2)+2
    nodo, tramo = nodos[inodo], tramos[itramo]
    if tramo.setj.bv_restri[1].get(): # uy prescrito
        valor= tramo.setj.valores[1]
        a_glob[iec] = valor
        prescindir[iec] = True
    if tramo.setj.bv_restri[2].get(): # fi prescrito
        valor= tramo.setj.valores[2]
        a_glob[iec+1] = valor
        prescindir[iec+1] = True

    
    # resolucion del sistema de ecuaciones
    # K_glob tiene 4*ntr x 4*ntr terminos, con los resortes incorporados
    # a_glob tiene 4*ntr terminos, con los valores de u no nulos puestos
    # f_glob tiene 4*ntr terminos, con las aportaciones de py Fext & Mext
    # prescindir[4*ntr] tiene True si el gdl esta restringido.
    # Uso auxiliares '_'
    
    # paso lo conocido al miembro dcho:
    f_= f_glob - np.matmul(K_glob, a_glob)
    # filtro f & K:
    quedar= np.logical_not(prescindir)
    f_= f_[quedar]
    K_= K_glob[quedar]
    K_= K_[:,quedar] # asi funciona, con K_glob[quedar, quedar] -> NO
                     # K_= K_glob[np.ix_(quedar, quedar)] -> SI funciona
    
    # resuelvo:
    try:
        a_= np.linalg.solve(K_, f_)
    except :
        texto1='Detectada configuración errónea.'
        texto2='''El sistema de ecuaciones del sub-problema de flexión no se puede resolver. Compruebe que no existen mecanismos locales en la configuración proporcionada.'''
        messagebox.showinfo(parent=v0, title='MAL', message=texto1, detail=texto2)
        return(1)
        
        
    # llevo el resultado a a_glob:
    i=0
    for igdl in range(4*ntramos):
        if quedar[igdl]:
            a_glob[igdl]=a_[i]
            i+=1
    # calculo reacciones en apoyos:
    f_apoyos = -np.matmul(K_glob, a_glob) + f_glob  # solo valen los de prescindir[]
    # las limpio para no inducir a error (ó las uso como comprob ~0?)
    #f_apoyos[quedar]=0

    for itramo in range (ntramos):

        # guardo en fa_tr[ [VMfia--]ppio, [VMfia--]fin  ] las fuerzas y desplazam

        tramo=tramos[itramo]
        tramo.fa_tr= [[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.]] 
        
        tramo.fa_tr[0][3] = a_glob[4*itramo]
        tramo.fa_tr[0][2] = a_glob[4*itramo+1] # el baile fi-u
        tramo.fa_tr[1][3] = a_glob[4*itramo+2] 
        tramo.fa_tr[1][2] = a_glob[4*itramo+3] 
        
        EI= tramo.ctes[1]*tramo.ctes[2]
        L=  tramo.ctes[0]
        a12, a6, a4, a2 = 12*EI/L**3, 6*EI/L**2, 4*EI/L, 2*EI/L
        K_tr= np.array([[a12, -a6, -a12, -a6],
                        [-a6,  a4,  a6,  a2],
                        [-a12, a6, a12, a6],
                        [-a6,  a2,  a6, a4]])
        p,q = tramo.carga[1][0],  tramo.carga[1][1]
        f0_tr=np.array([(-0.35*p-0.15*q)*L,
                        (p/20+q/30)*L*L,
                        (-0.15*p-0.35*q)*L,
                        (-p/30-q/20)*L*L     ])

        a_tr = a_glob[4*itramo:4*itramo+4]
        f_tr = np.matmul(K_tr, a_tr) + f0_tr 
            # es [Vi*, Mi*, Vj*, Mj*] del tramo (son fuerzas, no esf)
        tramo.fa_tr[0][0]=f_tr[0]
        tramo.fa_tr[0][1]=f_tr[1]
        tramo.fa_tr[1][0]=f_tr[2]
        tramo.fa_tr[1][1]=f_tr[3]


        # REACCIONES (AQUI ASIGNO LAS DE FLEXION)
        # guardo en tramo.reacc[  [ [],[],[] ], [ [],[],[] ]  ] iój / VMóN / vacio o valor
        # las reacc f_apoyos de apoyos & las de resortes (que calculo). 

        tramo.reacc=[   [ 0,0,0. ] , [ 0,0,0.]   ]
        igdl=4*itramo
        
        # reacc de apoyos:
        i=igdl
        if prescindir[i]: tramo.reacc[0][0]= f_apoyos[i]
        i=igdl+1
        if prescindir[i]: tramo.reacc[0][1]= f_apoyos[i]
        i=igdl+2
        if prescindir[i]: tramo.reacc[1][0]= f_apoyos[i]
        i=igdl+3
        if prescindir[i]: tramo.reacc[1][1]= f_apoyos[i]

        # reacc de resortes:
        if tramo.seti.bv_restri[3].get(): # ky
            tramo.reacc[0][0]= -tramo.seti.valores[3]*tramo.fa_tr[0][3]
        if tramo.setj.bv_restri[3].get():
            tramo.reacc[1][0]= -tramo.setj.valores[3]*tramo.fa_tr[1][3]
        if tramo.seti.bv_restri[4].get(): # kT
            tramo.reacc[0][1]= -tramo.seti.valores[4]*tramo.fa_tr[0][2]
        if tramo.setj.bv_restri[4].get():
            tramo.reacc[1][1]= -tramo.setj.valores[4]*tramo.fa_tr[1][2]



    # hacemos un calculo similar para los axiles
    
    # Hay 2 incog por tramo: (w)_i (w)_j
    K_globN= np.zeros((2*ntramos, 2*ntramos))
    a_globN= np.zeros(2*ntramos)
    f_globN= np.zeros(2*ntramos)
    prescindirN= np.array([False]*(2*ntramos))

    # Procedo por nodos, equilibrando cada uno. El primero y el ultimo van aparte.
    for inodo in range (1,nnodos-1):
        nodo=nodos[inodo]
        iec= 2*inodo-1
        itramo0, itramo1= inodo-1, inodo
        iin0, iin1 = 2*itramo0, 2*itramo1
        
        L0=tramos[itramo0].ctes[0]
        L1=tramos[itramo1].ctes[0]
        EAL0=tramos[itramo0].ctes[1]*tramos[itramo0].ctes[3]/L0
        EAL1=tramos[itramo1].ctes[1]*tramos[itramo1].ctes[3]/L1

        if nodo.bv_gdl[0].get():  # gdl 0 compartido
            K_globN[iec,iin0]  += -EAL0
            K_globN[iec,iin0+1]+=  EAL0
            K_globN[iec,iin1]  +=  EAL1
            K_globN[iec,iin1+1]+= -EAL1
            Fext= (tramos[itramo0].carga[2][1]+tramos[itramo1].carga[2][0] +
                   ( tramos[itramo0].carga[0][0]/6+tramos[itramo0].carga[0][1]/3)*L0 +
                   ( tramos[itramo1].carga[0][0]/3+tramos[itramo1].carga[0][1]/6)*L1 )
            f_globN[iec] += Fext 
            K_globN[iec+1,iin0+1]=  1.
            K_globN[iec+1,iin1]  = -1.

        else:                     # gdl 0 no compartido
            K_globN[iec,iin0]  += -EAL0
            K_globN[iec,iin0+1]+=  EAL0
            K_globN[iec+1,iin1]  +=  EAL1
            K_globN[iec+1,iin1+1]+= -EAL1
            Fext= (  tramos[itramo0].carga[2][1]+
                   ( tramos[itramo0].carga[0][0]/6+tramos[itramo0].carga[0][1]/3)*L0 )
            f_globN[iec]   += Fext
            Fext= (  tramos[itramo1].carga[2][0] +
                   ( tramos[itramo1].carga[0][0]/3+tramos[itramo1].carga[0][1]/6)*L1 )
            f_globN[iec+1] += Fext
            
    # El primer nodo:
    iec, iin = 0, 0
    EAL=tramos[0].ctes[1]*tramos[0].ctes[3]/tramos[0].ctes[0]
    L0=tramos[0].ctes[0]
    K_globN[iec,iin]  +=  EAL
    K_globN[iec,iin+1]+= -EAL
    f_globN[iec]   += ( tramos[0].carga[2][0] +
                        ( tramos[0].carga[0][0]/3+tramos[0].carga[0][1]/6)*L0 )

    # El ultimo nodo:
    iec, iin = 2*ntramos-1, 2*ntramos-2
    EAL=tramos[-1].ctes[1]*tramos[-1].ctes[3]/tramos[-1].ctes[0]
    L1 = tramos[-1].ctes[0]
    K_globN[iec,iin]  += -EAL
    K_globN[iec,iin+1]+=  EAL
    f_globN[iec]   += ( tramos[-1].carga[2][1] +
                        ( tramos[-1].carga[0][0]/6+tramos[-1].carga[0][1]/3)*L1 )

    # pongo ecuaciones prescindibles (& valor del gdl prescrito):

    for inodo in range(1, nnodos-1):
        nodo=nodos[inodo]
        tramo0, tramo1= tramos[inodo-1], tramos[inodo]
        iec= 2*inodo - 1

        if tramo0.setj.bv_restri[0].get(): # ux prescrito
            if nodo.bv_gdl[0].get(): # compartido
                valor= tramo0.setj.valores[0] or tramo1.seti.valores[0]
                a_globN[iec] = valor
                prescindirN[iec] = True
            else: # no compartido
                valor= tramo0.setj.valores[0]
                a_globN[iec] = valor
                prescindirN[iec] = True
                
        if tramo1.seti.bv_restri[0].get(): # ux prescrito
            if nodo.bv_gdl[0].get(): # esta compartido
                pass # dejar la ec de 1, -1 como esta (que calcule)
            else: # no esta compartido
                valor = tramo1.seti.valores[0]
                a_globN[iec+1] = valor
                prescindirN[iec+1] = True
        
    # el primer nodo:
    inodo, itramo, iec = 0, 0, 0
    nodo, tramo = nodos[inodo], tramos[itramo]
    if tramo.seti.bv_restri[0].get(): # ux prescrito
        valor= tramo.seti.valores[0]
        a_globN[iec] = valor
        prescindirN[iec] = True

    # el ultimo nodo:
    inodo, itramo, iec = nnodos-1, ntramos-1, 2*(ntramos)-1
    nodo, tramo = nodos[inodo], tramos[itramo]
    if tramo.setj.bv_restri[0].get(): # ux prescrito
        valor= tramo.setj.valores[0]
        a_globN[iec] = valor
        prescindirN[iec] = True


    # resolucion del sistema de ecuaciones
    # K_globN tiene 2*ntr x 2*ntr terminos. No admite resortes en x.
    # a_globN tiene 2*ntr terminos, con los valores de u no nulos puestos
    # f_globN tiene 2*ntr terminos, con las aportaciones de px & Fext
    # prescindirN[2*ntr] tiene True si el gdl esta restringido.
    # Uso auxiliares '_', las mismas que antes (se pierden, ok)
    
    # paso lo conocido al miembro dcho:
    f_= f_globN - np.matmul(K_globN, a_globN)
    # filtro f & K:
    quedarN= np.logical_not(prescindirN)
    f_= f_[quedarN]
    K_= K_globN[quedarN]
    K_= K_[:,quedarN] # asi funciona, con K_glob[quedar, quedar] -> NO
                      # K_= K_glob[np.ix_(quedar, quedar)] -> SI funciona
    
    # resuelvo:
    try:
        a_= np.linalg.solve(K_, f_)
    except :
        texto1='Detectada configuración errónea.'
        texto2='''El sistema de ecuaciones del sub-problema de tracc-compr. no se puede resolver. Compruebe que no existen mecanismos locales en la configuración proporcionada.'''
        messagebox.showinfo(parent=v0, title='MAL', message=texto1, detail=texto2)
        return(1)

    # llevo el resultado a a_glob:
    i=0
    for igdl in range(2*ntramos):
        if quedarN[igdl]:
            a_globN[igdl]=a_[i]
            i+=1
    # calculo reacciones en apoyos:
    f_apoyosN = np.matmul(K_globN, a_globN) - f_globN # solo valen los de prescindir[]
    # las limpio para no inducir a error (ó las uso como comprob ~0?)
    #f_apoyosN[quedarN]=0


    for itramo in range (ntramos):
        tramo=tramos[itramo]

        # guardo en fa_tr[ [x,x,x,x, N,w]ppio, [x,x,x,x, N,w]fin  ] las fuerzas y desplazam

        tramo.fa_tr[0][5] = a_globN[2*itramo]
        tramo.fa_tr[1][5] = a_globN[2*itramo+1] 
        
        EAL= tramo.ctes[1]*tramo.ctes[3]/tramo.ctes[0]

        K_tr= np.array([[1, -1], [-1, 1]]) * EAL
        r,s = tramo.carga[0][0],  tramo.carga[0][1]
        f0_tr=-np.array([ r/3+s/6 ,  r/6+s/3 ]) * tramo.ctes[0]

        a_tr = a_globN[2*itramo:2*itramo+2]
        f_tr = np.matmul(K_tr, a_tr) + f0_tr 
            # es [Fi*,Fj*] del tramo (son fuerzas, no esf)
        tramo.fa_tr[0][4]=f_tr[0]
        tramo.fa_tr[1][4]=f_tr[1]

        # REACCIONES (AQUI ASIGNO LAS DE AXIL)
        # guardo en tramo.reacc[  [ c,c,c ], [ c,c,c ]  ] iój / VMóN / vacio o valor
        # las reacc f_apoyosN de apoyos (no hay resortes N). 

        igdl=2*itramo
        i=igdl
        if prescindirN[i]: tramo.reacc[0][2]= f_apoyosN[i]
        i=igdl+1
        if prescindirN[i]: tramo.reacc[1][2]= f_apoyosN[i]


    # Siguen varias asignaciones para ambos (flexion & axil): 
    # - funciones poly en los tramos
    # - ceros & maximos en los tramos

    # asigno las funciones en los tramos:
    for tramo in tramos: tramo.ponfunciones()

    # localizar ceros de V,M,fi,u,N,w interiores al tramo, si hay. Deja en tr.haycero[]
    # seis sublistas, bien vacias o bien con los x de los ceros:

    for itramo in range(ntramos):
        tramo=tramos[itramo]
        for igr in range(6): # para las 6 graficas
            tramo.haycero[igr]=[]
            funcion=tramo.funciones[igr] # en forma de polinomios
            ceros= funcion.roots()
            if np.isreal(any(ceros)):
                for c in ceros:
                    if np.isreal(c):
                        if 0.04*tramo.ctes[0] < c < 0.96*tramo.ctes[0]:
                            tramo.haycero[igr].append(np.real(c))

    # localizar maximos de V,M,fi,u,N,w interiores al tramo, si hay. Deja en tr.maximos[]
    # seis sublistas, bien vacias o bien con duplas tipo [x_max,valor],[x_max,valor]... :

    for itramo in range(ntramos):
        tramo=tramos[itramo]
        tramo.maximos= [ [],[],[],[],[],[]  ]

        for igr in range(6): # para las 6 graficas

            haymax, xdemax = False, []
            if igr==0: # el del cortante hay que calcularlo ahora:
                if ((tramo.carga[1][0]>0. and tramo.carga[1][1]<0) or
                    (tramo.carga[1][0]<0. and tramo.carga[1][1]>0)):
                    haymax=True
                    xdemax.append(  abs(tramo.carga[1][0])*tramo.ctes[0]/
                              (abs(tramo.carga[1][0])+abs(tramo.carga[1][1]) )  )
            elif igr==4: # el de N tambien
                if ((tramo.carga[0][0]>0. and tramo.carga[0][1]<0) or
                    (tramo.carga[0][0]<0. and tramo.carga[0][1]>0)):
                    haymax=True
                    xdemax.append(  abs(tramo.carga[0][0])*tramo.ctes[0]/
                              (abs(tramo.carga[0][0])+abs(tramo.carga[0][1]) )  )
            else:
                haymax=bool(tramo.haycero[igr-1])
                xdemax=tramo.haycero[igr-1] # puede haber varios, es una lista

            if haymax:
                for xm in xdemax:
                    if not(np.isclose(xm, tramo.ctes[0])) and xm>tramo.ctes[0]/1.e3:
                        ym= tramo.funciones[igr](xm)
                        tramo.maximos[igr].append([xm,ym])

    salida_texto()
    salida_grafica()

    return
    

    
######### EMPIEZA EL PROGRAMA  ###########



hard_ntramos=9
v3,v4, frame_lienzoN, frame_lienzoM = None, None, None, None
tramos, nodos, ntramos, nnodos, scal_px, ltotal = [], [], 0, 0, 0, 0
btn_do_cargas, btn_do_apoyos, btn_do_nonulos= '', '', ''
entries_vistas=False

# esto es para que encuentre 'arte' cuando se corre desde terminal / lanzador
try:
    os.chdir(os.path.dirname(__file__)) 
    print('\nEl directorio de trabajo es:\n',os.getcwd())
except FileNotFoundError:
    pass


###########################
### INTERFAZ DE USUARIO ###
###########################



v0=Tk()
v0.title('JContinua')

estilo = ttk.Style()
estilo.configure('jc_blue.TCheckbutton', background='#D7ECFF')
estilo.configure('jc_green.TCheckbutton', background='#DAFFD7')
estilo.configure('jc.TLabelframe.Label', foreground ='green')
estilo.configure('jc_blue.TFrame', background='#D7ECFF')  # no lo uso?
estilo.configure('jc_white.TLabel', background='white')
estilo.configure('jc_grey.TEntry', background='#EDECEB')

# figs aqui pq hay bug en Tk, "mantener referencia a objetos entre llamadas"
figs={
'nuxuyfi': PhotoImage( file= './arte/nuxuyfi.png'),
'nuy': PhotoImage(     file= './arte/nuy.png'),
'nuxfi': PhotoImage(   file= './arte/nuxfi.png'),
'nuxuy': PhotoImage(   file= './arte/nuxuy.png'),
'nux': PhotoImage(     file= './arte/nux.png'),
'nuyfi': PhotoImage(   file= './arte/nuyfi.png'),
'suxd': PhotoImage(    file= './arte/suxd.png'),
'suxi': PhotoImage(    file= './arte/suxi.png'),
'suxuy': PhotoImage(   file= './arte/suxuy.png'),
'skt': PhotoImage(     file= './arte/skt.png'),
'suxfid': PhotoImage(  file= './arte/suxfid.png'),
'suxuyfid': PhotoImage(file= './arte/suxuyfid.png'),
'suy': PhotoImage(     file= './arte/suy.png'),
'sfi': PhotoImage(     file= './arte/sfi.png'),
'sky': PhotoImage(     file= './arte/sky.png'),
'suxfii': PhotoImage(  file= './arte/suxfii.png'),
'suxuyfii': PhotoImage(file= './arte/suxuyfii.png'),
'suyfi':PhotoImage(    file= './arte/suyfi.png'),
'dgn':PhotoImage(    file= './arte/dgn.png'),
'dgp':PhotoImage(    file= './arte/dgp.png'),
'dmn':PhotoImage(    file= './arte/dmn.png'),
'dmp':PhotoImage(    file= './arte/dmp.png'),
'dvn':PhotoImage(    file= './arte/dvn.png'),
'dvp':PhotoImage(    file= './arte/dvp.png'),
'dnp':PhotoImage(    file= './arte/dnp.png'),
'dnn':PhotoImage(    file= './arte/dnn.png'),
'dwp':PhotoImage(    file= './arte/dwp.png'),
'dwn':PhotoImage(    file= './arte/dwn.png')  }

frame_general = ttk.Frame(v0, padding='4')
frame_general.grid(sticky=(N, W, E, S))

frame_botones= ttk.Frame(frame_general, width=1100,height=50)
frame_botones.grid(column=0,row=0,ipadx=3)

frame_lienzo= Canvas(frame_general, width=1100, height=550, background='white')
frame_lienzo.grid(column=0, row=1, sticky=(N, W, E, S))


# rellenar frame_botones (tira superior en la ventana):
btn_tramos= ttk.Button(frame_botones, text='tramos', command=di_tramos)
btn_tramos.grid(column=0, row=0, padx=3, pady=5, sticky='w')

btn_calcular= ttk.Button(frame_botones, text='calcular', command=calcular)
btn_calcular.grid(column=1, row=0, padx=3, pady=5, sticky='w')

btn_imprimir= ttk.Button(frame_botones, text='imprimir', command=imprimir)
btn_imprimir.grid(column=2, row=0, padx=3, pady=5, sticky='w')

btn_ayuda= ttk.Button(frame_botones, text='ayuda', command=ayuda)
btn_ayuda.grid(column=3, row=0, padx=3, pady=5, sticky='e')

btn_salir= ttk.Button(frame_botones, text='salir', command=salir)
btn_salir.grid(column=4, row=0, padx=3, pady=5, sticky='e')

v0.protocol('WM_DELETE_WINDOW', salir)
v0.mainloop()














