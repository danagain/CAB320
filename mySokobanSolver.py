'''

The partially defined functions and classes of this module 
will be called by a marker script. 

You should complete the functions and classes according to their specified interfaces.


'''

import search

import sokoban
from sokoban import Warehouse


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def my_team():
    '''
    Return the list of the team members of this assignment submission as a list
    of triplet of the form (student_number, first_name, last_name)

    '''
    return [(9493671, 'Daniel', 'Huffer'), (1234568, 'Grace', 'Hopper'), (1234569, 'Eva', 'Tardos')]
    raise NotImplementedError()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def taboo_cells(warehouse):
    '''
    Identify the taboo cells of a warehouse. A cell is called 'taboo'
    if whenever a box get pushed on such a cell then the puzzle becomes unsolvable.
    When determining the taboo cells, you must ignore all the existing boxes,
    simply consider the walls and the target  cells.
    Use only the following two rules to determine the taboo cells;
     Rule 1: if a cell is a corner and not a target, then it is a taboo cell.
     Rule 2: all the cells between two corners along a wall are taboo if none of
             these cells is a target.

    @param warehouse: a Warehouse object

    @return
       A string representing the puzzle with only the wall cells marked with
       an '#' and the taboo cells marked with an 'X'.
       The returned string should NOT have marks for the worker, the targets,
       and the boxes.
    '''
    ##         "INSERT YOUR CODE HERE"


    #problem_file = "./warehouses/warehouse_03.txt"  ##Load up maps just for testing
    #warehouse.read_warehouse_file(problem_file)

    ######This code is from the Sokoban file to get coords of objects in the warehouse


    freetile = list()
    X, Y = zip(*warehouse.walls)
    x_size, y_size = 1 + max(X), 1 + max(Y)
    vis = [[" "] * x_size for y in range(y_size)]
    for (x, y) in warehouse.walls:
        vis[y][x] = "#"
    for (x, y) in warehouse.targets:
        vis[y][x] = "."
    #######
    stringWarehouse = "\n".join(["".join(line) for line in vis])
    print(stringWarehouse)
    colx = len(vis[1])  # Length of map in columns  (x)
    rowy = len(vis) - 1  # Length of map in rows (y)
    tabooTiles = list()  # List for bad tiles

    # scan entire warehouse taking coords for every blank tile
    for x in range(colx):
        for y in range(rowy):
            if vis[y][x] != "#" or vis[y][x] != ".":
                freetile.append((x, y))



                # Take all the corner tiles ( Takes blank spots outside the map, not sure if an issue)
    for tiles in freetile:
        tile_x, tile_y = tiles  # assign coords to all the blank tiles
        north_x, north_y = tile_x, tile_y - 1  # set up north facing position
        south_x, south_y = tile_x, tile_y + 1  # south ..
        east_x, east_y = tile_x + 1, tile_y  # east..
        west_x, west_y = tile_x - 1, tile_y  # west..
        if ((north_x, north_y) in warehouse.walls and (east_x, east_y) in warehouse.walls) \
                or \
                ((south_x, south_y) in warehouse.walls and (east_x, east_y) in warehouse.walls) \
                or \
                ((south_x, south_y) in warehouse.walls and (west_x, west_y) in warehouse.walls) \
                or \
                ((north_x, north_y) in warehouse.walls and (west_x, west_y) in warehouse.walls):
            tabooTiles.append(tiles)

    for i in range(len(tabooTiles) - 1):
        t_x, t_y = tabooTiles[i]
        if vis[t_y][t_x] != ".":
            vis[t_y][t_x] = "X"

            #### This code should cover all bottom walls between to corners ####


            # For the length of the array of " " spaces
            for j in range(len(tabooTiles) - 1):
                j_x, j_y = tabooTiles[j]  # Take the coordinates
                for i in range(len(tabooTiles) - 2):  # for the length of tabooTiles -1 start one element ahead
                    i_x, i_y = tabooTiles[i + 1]
                    if i_y == j_y and i_x != j_x:  # if y coords equal at some point during both loops
                        coverwall = None
                        counter = 0  # and the element infront is a hash tag
                        for k in range(j_x, i_x):  # for the range of the x coords of the matching y coords
                            if vis[i_y - 1][k] == "#" or vis[i_y + 1][k] == "#":  # loop across the x range checking it's a wall
                                counter = counter + 1  # count number of blocks
                                if counter == abs(j_x - i_x):  # if the counter = the range of the x coords
                                    for p in range(j_x, i_x):  # then loop across marking tabooTiles
                                        if vis[j_y][p] == ".":
                                            coverwall = False
                        if coverwall == True:
                            for p2 in range(j_x, i_x):  # then loop across marking tabooTiles
                                vis[j_y][p2] = "X"



                                            #### This code should cover all side walls between to corners ####

                                        # For the length of the array of " " spaces
        for j in range(len(tabooTiles) - 1):
            j_x, j_y = tabooTiles[j]  # Take the coordinates
            for i in range(len(tabooTiles) - 2):  # for the length of tabooTiles -1 start one element ahead
                i_x, i_y = tabooTiles[i + 1]
                if i_x + 1 < colx:
                    if i_x == j_x and i_y != j_y:  # if y coords equal at some point during both loops
                        coverwall = None
                        counter = 0  # and the element infront is a hash tag
                        for k in range(j_y, i_y):  # for the range of the x coords of the matching y coords
                            if (vis[k][i_x - 1] == "#" or vis[k][
                                    i_x + 1] == "#"):  # loop across the x range checking it's a wall
                                counter = counter + 1  # count number of blocks
                                if counter == abs(j_y - i_y):  # if the counter = the range of the x coords
                                    for p in range(j_y, i_y):  # then loop across marking tabooTiles
                                        if vis[p][j_x] == ".":
                                            coverwall = False
                                        else: coverwall = True
                        if coverwall == True:
                            for p2 in range(j_y, i_y):
                                vis[p2][j_x] = "X"
    #wall_coords = []
    #for (x, y) in warehouse.walls:
    #   wall_coords.append((x,y))
    #print( wall_coords)

    ##CLear up "X" outside of map (almost works perfectly except a few maps)
    for y in range(rowy):
        Xdata = []
        Wdata = []
        for x in range(colx):
            if vis[y][x] == "X":
                Xdata.append(x)
            if vis[y][x] == "#":
                Wdata.append(x)
        Xdata.sort()
        Wdata.sort()
        if len(Xdata) >= 1 and len(Wdata) >= 1:
            if max(Xdata) > max(Wdata):
                vis[y][max(Xdata)] = " "
        if len(Xdata) >= 1 and len(Wdata) >= 1:
            if min(Xdata) < min(Wdata):
                vis[y][min(Xdata)] = " "

    for (x, y) in warehouse.walls:
        vis[y][x] = "#"
    for (x, y) in warehouse.targets:
        vis[y][x] = " "
    stringWarehouse = "\n".join(["".join(line) for line in vis])
    print(stringWarehouse)
    return stringWarehouse


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class SokobanPuzzle(search.Problem):
    '''
    Class to represent a Sokoban puzzle.
    Your implementation should be compatible with the
    search functions of the provided module 'search.py'.

    	Use the sliding puzzle and the pancake puzzle for inspiration!

    '''

    ##         "INSERT YOUR CODE HERE"

    def __init__(self,warehouse):
        self.initial_state = self.LoadProblem(warehouse)
        self.current_state = self.initial_state

    def LoadProblem(self, filePath):
        environment = Warehouse()
        environment.read_warehouse_file(filePath)
        return environment

    def getState(self, state=None):
        if state == None:
            state = self.current_state
        return state


    def actions(self, state):
        """
        Return the list of actions that can be executed in the given state
        if these actions do not push a box in a taboo cell.
        The actions must belong to the list ['Left', 'Down', 'Right', 'Up']
        """
        raise NotImplementedError


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def check_action_seq(warehouse, action_seq):
    '''

    Determine if the sequence of actions listed in 'action_seq' is legal or not.

    Important notes:
      - a legal sequence of actions does not necessarily solve the puzzle.
      - an action is legal even if it pushes a box onto a taboo cell.

    @param warehouse: a valid Warehouse object

    @param action_seq: a sequence of legal actions.
           For example, ['Left', 'Down', Down','Right', 'Up', 'Down']

    @return
        The string 'Failure', if one of the action was not successul.
           For example, if the agent tries to push two boxes at the same time,
                        or push one box into a wall.
        Otherwise, if all actions were successful, return
               A string representing the state of the puzzle after applying
               the sequence of actions.  This must be the same string as the
               string returned by the method  Warehouse.__str__()
    '''

    ##         "INSERT YOUR CODE HERE"

    raise NotImplementedError()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def solve_sokoban_elem(warehouse):
    '''
    This function should solve using elementary actions
    the puzzle defined in a file.

    @param warehouse: a valid Warehouse object

    @return
        A list of strings.
        If puzzle cannot be solved return ['Impossible']
        If a solution was found, return a list of elementary actions that solves
            the given puzzle coded with 'Left', 'Right', 'Up', 'Down'
            For example, ['Left', 'Down', Down','Right', 'Up', 'Down']
            If the puzzle is already in a goal state, simply return []
    '''

    ##         "INSERT YOUR CODE HERE"

    raise NotImplementedError()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def can_go_there(warehouse, dst):
    '''
    Determine whether the worker can walk to the cell dst=(row,col)
    without pushing any box.

    @param warehouse: a valid Warehouse object

    @return
      True if the worker can walk to cell dst=(row,col) without pushing any box
      False otherwise
    '''

    ##         "INSERT YOUR CODE HERE"

    raise NotImplementedError()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def solve_sokoban_macro(warehouse):
    '''
    Solve using macro actions the puzzle defined in the warehouse passed as
    a parameter. A sequence of macro actions should be
    represented by a list M of the form
            [ ((r1,c1), a1), ((r2,c2), a2), ..., ((rn,cn), an) ]
    For example M = [ ((3,4),'Left') , ((5,2),'Up'), ((12,4),'Down') ]
    means that the worker first goes the box at row 3 and column 4 and pushes it left,
    then goes the box at row 5 and column 2 and pushes it up, and finally
    goes the box at row 12 and column 4 and pushes it down.

    @param warehouse: a valid Warehouse object

    @return
        If puzzle cannot be solved return ['Impossible']
        Otherwise return M a sequence of macro actions that solves the puzzle.
        If the puzzle is already in a goal state, simply return []
    '''

    ##         "INSERT YOUR CODE HERE"

    raise NotImplementedError()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if __name__ == "__main__":
    a = SokobanPuzzle("./warehouses/warehouse_03.txt")
    b = a.getState()
    taboo_cells(b)
    #print(b)
    #taboo_cells(wh)
