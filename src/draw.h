#ifndef DRAW_H
#define DRAW_H

typedef struct __Tdraw {
  int             *xcoord; /* node coordinates on the x axis */
  int             *ycoord; /* node coordinates on the y axis */
  int          page_width;
  int         page_height;
  int      tree_box_width;

  double max_dist_to_root;
}tdraw;



void Dist_To_Root_Pre(node *a, node *d, edge *b, arbre *tree);
void Dist_To_Root(edge *b_root, arbre *tree);
void Get_X_Coord_Pre(node *a, node *d, edge *b, tdraw *w, arbre *tree);
void Get_X_Coord(edge *b_root, tdraw *w, arbre *tree);
tdraw *Make_Tdraw_Struct(arbre *tree);
void Init_Tdraw_Struct(tdraw *d);
void Get_Tree_Box_Width(tdraw *w, arbre *tree);
void Get_Y_Coord_Post(node *a, node *d, edge *b, int *next_y_slot, tdraw *w, arbre *tree);
void Get_Y_Coord(edge *b_root, tdraw *w, arbre *tree);
void Draw_Tree(edge *b_root, tdraw *w, arbre *tree);
double Get_Max_Dist_To_Root(arbre *tree);
void Print_Tree_Postscript(edge *b_root, FILE *fp, int tree_num, tdraw *w, arbre *tree);
void Print_Tree_Postscript_Pre(node *a, node *d, FILE *fp, tdraw *w, arbre *tree);
void Print_Postscript_EOF(FILE *fp);
void Print_Postscript_Header(int n_pages, FILE *fp);


#endif
