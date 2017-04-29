
'''

The partially defined functions and classes of this module
will be called by a marker script.

You should complete the functions and classes according to their specified interfaces.


'''

import search
from search import depth_first_graph_search
from search import breadth_first_graph_search
import itertools
import operator
from sokoban import Warehouse

actions_to_goal = []
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def my_team():
    '''
    Return the list of the team members of this assignment submission as a list
    of triplet of the form (student_number, first_name, last_name)

    '''
    return [ (9448551, 'Daniel', 'Huffer'), (1234568, 'Grace', 'Hopper'), (1234569, 'Eva', 'Tardos') ]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


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

    freetile = list() #making a list to store the coordinates of tiles with nothing on them
    X, Y = zip(*warehouse.walls) #finding all the coordinates of the walls
    x_size, y_size = 1 + max(X), 1 + max(Y)
    vis = [[" "] * x_size for y in range(y_size)]
    for (x, y) in warehouse.walls:
        vis[y][x] = "#"
    for (x, y) in warehouse.targets:
        vis[y][x] = "."
    colx = len(vis[1])  # Length of map in columns  (x)
    rowy = len(vis) - 1  # Length of map in rows
    taboo_tiles = list()  # List for bad tiles
    # scan entire warehouse taking coords for every blank tile
    for x in range(colx):
        for y in range(rowy):
            if vis[y][x] != "#" or vis[y][x] != ".":
                freetile.append((x, y))
    # Take all the corner tiles ( Takes blank spots outside the map needs to be resolved)
    corner_taboos(freetile, warehouse, taboo_tiles, vis)
    taboo_coordinates = []
    for i in range(len(taboo_tiles) - 1):
        t_x, t_y = taboo_tiles[i]
        if vis[t_y][t_x] != ".":
            vis[t_y][t_x] = "X"
    taboo_walls(taboo_tiles, vis , colx)
    ## Need to find a way to find all the clear tiles inside the warehouse
    ## Tutor said to do a breadth_first_graph_search with no goal state and save all the positions

    ## Taboo_tiles is working except for the fact we need to remove some of the taboo_tiles that
    ## are popping up outside the map.

    stringFormat = clear_outter_taboos(rowy, colx, vis, warehouse)
    return stringFormat

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class SokobanPuzzle(search.Problem):


    '''
    Class to represent a Sokoban puzzle.
    Your implementation should be compatible with the
    search functions of the provided module 'search.py'.

    	Use the sliding puzzle and the pancake puzzle for inspiration!

    '''
    
    def __init__(self, warehouse):
        self.warehouse = warehouse # assign the warehouse string to instance var warehouse
        self.initial = (warehouse.worker,tuple(warehouse.boxes)) # init state is dbl tuple of worker and boxes
        self.taboos = taboo_coordinates(warehouse) # load up all taboo coordinates to instance variable taboos

    """
    goal_test is checking all the target coordinates match the set of box coordinates
    """
    def goal_test(self, state):
        return set(self.warehouse.targets) == set(state[1])  # all boxes match target locations

    def actions(self, state):
        """
        Return the list of actions that can be executed in the given state
        if these actions do not push a box in a taboo cell.
        The actions must belong to the list ['Left', 'Down', 'Right', 'Up']
        """
        worker = state[0]
        boxes = state[1]
        walls = self.warehouse.walls
        taboo = self.taboos
        valid_move = []
        moves = {"Up":(worker[0],worker[1]-1), "Down": (worker[0],worker[1]+1), "Left": (worker[0] -1, worker[1]), "Right": (worker[0] + 1, worker[1])}  # getting dictionary errors on this line ?!??!
        for direction in moves:
            pos = moves[direction]
            if pos not in walls and pos not in boxes:
                valid_move.append(str(direction))
            if pos in boxes:
                if direction == "Up":
                    pos = (pos[0], pos[1] -1)
                if direction == "Down":
                    pos = (pos[0], pos[1] +1)
                if direction == "Left":
                    pos = (pos[0]-1, pos[1])
                if direction == "Right":
                    pos = (pos[0] + 1, pos[1])

                    if pos not in list(boxes)+taboo+walls:
                        valid_move.append(str(direction))
        return valid_move

    def print_solution(self, goal_node):
        if goal_node != None:
            path = goal_node.path()
            return self.goal_path(path)

    def goal_path(self, path):
        actions = []
        for node in path:
            if node.action is not None:
                actions.append(node.action)
        return actions


    # Results of the action
    def result(self, state, action):
        """
        Return theh state that results from executing the given
        action in the given state. The action must be one of
        self.actions(state).
        """
        #assert action in self.actions(state)
        worker = state[0]
        boxes = state[1]
        box_list = list(boxes)
        Actions = [action]
        for direction in Actions:
            #worker =  moves.get(direction) #  Should add the tuples together and assign to worker
            if direction == "Up":
                worker = (worker[0], worker[1]-1)
                actions_to_goal.append("Up")
            if direction == "Down":
                worker = (worker[0], worker[1]+1)
                actions_to_goal.append("Down")
            if direction == "Left":
                worker = (worker[0]-1, worker[1])
                actions_to_goal.append("Left")
            if direction == "Right":
                worker = (worker[0]+1, worker[1])
                actions_to_goal.append("Right")
            #  Player should move and if he moves into a box the box then needs to move a square in the same direction
            for box in boxes:
                b = box
                if box == worker:
                    if direction == "Up":
                        b = (b[0] , b[1]-1)
                    if direction == "Down":
                        b = (b[0], b[1]+1)
                    if direction == "Left":
                        b = (b[0]-1, b[1])
                    if direction == "Right":
                        b = (b[0]+1, b[1])
                            # move the box assign to a temp variable
                    box_list.remove(box)  # remove the current location of the box
                    box_list.append(b)   # add the location of the new box position
            continue  #  break the loop as the box is moved
        return (worker,tuple(box_list))

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
    walls = warehouse.walls
    worker = warehouse.worker  # Workers location
    box = warehouse.boxes  # List tuples containing box coordinates
    directions = ["Up", "Down", "Left", "Right"]
    move = {"Up": (worker[0], worker[1]-1), "Down": (worker[0], worker[1]+1), \
            "Left": (worker[0]-1, worker[1]), "Right": (worker[0]+1, worker[1])}
    moves = {"Up": (0, -1), "Down": (0, 1), "Left": (-1, 0), "Right": (1, 0)}
    move_box = {"Up": (worker[0],worker[1]-2),"Down":(worker[0],worker[1]+2), \
                "Left": (worker[0]-2,worker[1]),"Right": (worker[0] + 2, worker[1])}
    for action in action_seq:
        for d in directions:
            if action == d:
                if (move[d] in walls) or (move[d] in box and move_box[d] in walls) or \
                        (move[d] in box and move_box[d] in box):
                    return 'Failure'
                else:
                    for b in box:
                        if b == move[d]:
                            b = move_box_func[d]
                    worker = move_player(worker, moves[d])
    string_warehouse = warehouse.copy(worker, box)
    return string_warehouse.__str__()


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
    puzzle = SokobanPuzzle(warehouse)
    sol_ts = breadth_first_graph_search(puzzle)
    #sol_ts = depth_first_graph_search(puzzle)
    ans = puzzle.print_solution(sol_ts)
    return ans


    #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""
Logic:
Sif the worker is on the y axis
for the range of the movement to the destination in the y axis, if the x axis meets a box
then return false, likewise for moving in the x direction return false if box is on same y coordinate
"""
def can_go_there(warehouse, dst):
    '''
    Determine whether the worker can walk to the cell dst=(row,col)
    without pushing any box.

    @param warehouse: a valid Warehouse object

    @return
      True if the worker can walk to cell dst=(row,col) without pushing any box
      False otherwise
    '''

    X, Y = zip(*warehouse.walls)
    x_size, y_size = 1 + max(X), 1 + max(Y)
    x_movement = []
    for x in range(warehouse.worker[0],dst[0]):
        x_movement.append(x)
    y_movement = []
    for y in range(warehouse.worker[1], dst[1]):
        y_movement.append(y)

    for b in warehouse.boxes: # Check every box
        if b[0] in x_movement and b[1] == warehouse.worker[1]: # If the box has the same x as the box, check y matches
            return False
        if b[1] in y_movement and b[0] == warehouse.worker[0]: # If the box y is in the y list and the x matches
            return False
        if dst[0] > x_size or dst[0] < 0: # If the destination is off the map
            return False
        if dst[1] > y_size or dst[1] < 0: # If the destination is off the map
            return False
        else:
            return True


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

    # use the actions defined for elementary actions but change all actions that do not push a block
    # into the coordinate that would result from those actions

    raise NotImplementedError()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"""
This function uses the itertools to find two matching y coordinates with taboo tiles
If there x values vary and there is no target points between them it will fill across the wall with taboo tiles
same process is repeated for matching x coordinates
"""
def taboo_walls(taboo_tiles, vis, colx):
    def check_top_bot_walls():
            wall_t = 0
            wall_b = 0
            for x in range(x1+1, x2):
                if vis[y1-1][x] == "#" and vis[y1][x] != ".":
                    wall_t += 1
                if vis[y1+1][x] == "#" and vis[y1][x] != ".":
                    wall_b += 1
                if wall_t == abs(x2 - (x1+1)):
                        for wall in range((x1+1), x2):
                            vis[y1][wall] = "X"
                if wall_b == abs(x2 - (x1+1)):
                        for w in range((x1+1), x2):
                            vis[y1][w] = "X"
    """
    The side walls is the same deal with bot and top walls function above, when it finds two matching x coordinates
    it will tile the wall with taboos if the y points vary and there is a confirmed wall
    """
    def check_side_walls():
        wall_l = 0
        wall_r = 0
        for y in range(y1+1, y2):
            if vis[y][x1 - 1] == "#" and vis[y][x1] != ".":
                wall_l += 1
            if wall_l == abs(y2 - (y1+1)):
                for w in range((y1 + 1), y2):
                   vis[w][x1] = "X"
            if x1 + 1 < colx:
                if vis[y][x1 + 1] == "#" and vis[y][x1] != ".":
                    wall_r += 1
                if wall_r == abs(y2 - (y1+1)):
                    for w in range((y1+1), y2):
                        vis[w][x1] = "X"
    for t1, t2 in itertools.permutations(taboo_tiles, 2):
        x1,y1 = t1
        x2,y2 = t2
        if y1 == y2:
            check_top_bot_walls()
        if x1 == x2:
            check_side_walls()
"""
This function checks every empty space looking to see if there is a corner (imagine standing staring into a right angle)
"""
def corner_taboos(freetile, warehouse, taboo_tiles, vis):
    for tiles in freetile: #Check every single free tile against every possible angle a wall could be
        tile_x, tile_y = tiles  # assign coords to all the blank tiles
        up_x, up_y = tile_x, tile_y - 1  # set up up facing position
        down_x, down_y = tile_x, tile_y + 1  # down ..
        right_x, right_y = tile_x + 1, tile_y  # right..
        left_x, left_y = tile_x - 1, tile_y  # left..
        if ((up_x, up_y) in warehouse.walls and (right_x, right_y) in warehouse.walls) and vis[tile_y][tile_x] != "." \
                or ((down_x, down_y) in warehouse.walls and (right_x, right_y) in warehouse.walls) \
                and vis[tile_y][tile_x] != "." or ((down_x, down_y) in warehouse.walls and (left_x, left_y) in\
                warehouse.walls) and vis[tile_y][tile_x] != "." or ((up_x, up_y) in warehouse.walls and
                (left_x, left_y) in warehouse.walls) and vis[tile_y][tile_x] != ".":
            taboo_tiles.append(tiles)

"""
## This function doesnt remove all the outter taboos, it is impossible to do so
# without finding all the free tiles strictly inside the walls first through a breadth first search or something first
"""
def clear_outter_taboos(rowy, colx, vis, warehouse):
    # Clear up "X" outside of map (almost works perfectly except a few maps)
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
    return stringWarehouse

"""
# This function is basically the same as the taboo_cells except it returns the coordinates
# of all the taboo cells instead of a multi-line string
"""
def taboo_coordinates(warehouse):
    freetile = list()  # making a list to store the coordinates of tiles with nothing on them
    X, Y = zip(*warehouse.walls)  # finding all the coordinates of the walls
    x_size, y_size = 1 + max(X), 1 + max(Y)
    vis = [[" "] * x_size for y in range(y_size)]
    for (x, y) in warehouse.walls:
        vis[y][x] = "#"
    for (x, y) in warehouse.targets:
        vis[y][x] = "."
    colx = len(vis[1])  # Length of map in columns  (x)
    rowy = len(vis) - 1  # Length of map in rows 󰀀
    taboo_tiles = list()  # List for bad tiles
    # scan entire warehouse taking coords for every blank tile
    for x in range(colx):
        for y in range(rowy):
            if vis[y][x] != "#" or vis[y][x] != ".":
                freetile.append((x, y))
    # Take all the corner tiles ( Takes blank spots outside the map needs to be resolved)
    corner_taboos(freetile, warehouse, taboo_tiles, vis)
    taboo_coords = []
    for i in range(len(taboo_tiles) - 1):
        t_x, t_y = taboo_tiles[i]
        if vis[t_y][t_x] != ".":
            vis[t_y][t_x] = "X"
    taboo_walls(taboo_tiles, vis, colx)
    clear_outter_taboos(rowy, colx, vis, warehouse)
    for x in range(colx):
        for y in range(rowy):
            if vis[y][x] == "X":
                taboo_coords.append((x, y))
    return taboo_coords


"""
simple functions to move the player and the box
by the addition of two tuples a and b
"""
def move_player(a,b):
    new = tuple(map(operator.add, a, b))
    return new

def move_box_func(a,b):
    new = tuple(map(operator.add, a, b))
    return new


"""
Try and find a way to test the search algorithm on
multiple warehouses  -- not working at the moment
"""
def test_puzzle(warehouse):
    puzzle = SokobanPuzzle(warehouse)
    ans = search.breadth_first_graph_search(puzzle)
    print(ans)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if __name__ == "__main__":
    print("main")


