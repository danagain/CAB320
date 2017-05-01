
'''

The partially defined functions and classes of this module
will be called by a marker script.

You should complete the functions and classes according to their specified interfaces.


'''

import search
from search import breadth_first_graph_search
import itertools
import operator
from sokoban import Warehouse
import time

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def my_team():
    '''
    Return the list of the team members of this assignment submission as a list
    of triplet of the form (student_number, first_name, last_name)

    '''
    return [ (9448551, 'Daniel', 'Huffer'), (7386044, 'Jenyfer', 'Florentina'), (9469184, 'Mitchell', 'Hurst') ]

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
    corner_taboos(freetile, warehouse, taboo_tiles, vis)  # Apply all corner taboos
    for i in range(len(taboo_tiles) - 1):
        t_x, t_y = taboo_tiles[i]
        if vis[t_y][t_x] != ".":
            vis[t_y][t_x] = "X"
    taboo_walls(taboo_tiles, vis , colx)  # Apply all wall taboos
    stringFormat = clear_outter_taboos(rowy, colx, vis, warehouse)  #  clean up warehouse
    return stringFormat

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


class SokobanPuzzle(search.Problem):
    '''
    Class to represent a Sokoban puzzle.
    Compatible with the search functions located in the
    'search.py' module
    '''

    def __init__(self, warehouse):
        """
        Class constructor for initialisation
        @param state: Warehouse object
        """
        self.warehouse = warehouse # assign the warehouse string to instance var warehouse
        self.initial = (warehouse.worker,tuple(warehouse.boxes)) # init state is dbl tuple of worker and boxes
        self.taboos = taboo_coordinates(warehouse) # load up all taboo coordinates to instance variable taboos

    def goal_test(self, state):
        """
        Determines if the current warehouse problem is in a goal state
        @param state: Current state of the warehouse

        @return
        True if all unique box coordinates match unique target coordinates
        Else, False
        """
        return set(self.warehouse.targets) == set(state[1])  # all boxes match target locations

    def actions(self, state):
        """
        Determines a list of viable worker actions based on current warehouse state

        @param state: Current state of the warehouse

        @return
        List of possible actions from current position e.g  ["Up", "Left"]
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
                if direction == "Up" and (pos[0],pos[1]-1) not in boxes and (pos[0],pos[1]-1) not in taboo and (pos[0],pos[1]-1) not in walls:
                    valid_move.append(str(direction))
                if direction == "Down"and (pos[0],pos[1]+1) not in boxes and (pos[0],pos[1]+1) not in taboo and (pos[0],pos[1]+1) not in walls:
                    valid_move.append(str(direction))
                if direction == "Left" and (pos[0]-1,pos[1]) not in boxes and (pos[0]-1,pos[1]) not in taboo and (pos[0]-1,pos[1]) not in walls:
                    valid_move.append(str(direction))
                if direction == "Right" and (pos[0]+1,pos[1]) not in boxes and (pos[0]+1,pos[1]) not in taboo and (pos[0]+1,pos[1]) not in walls:
                    valid_move.append(str(direction))
        return valid_move

    def print_solution(self, goal_node):
        """
        Print the sequence of actions resulting in the goal state
        if goal state was found.
        @param goal_node: The final node popped off the frontier in the search
        algorithm

        @return
        A string list of each action leading to the goal state
        """
        if goal_node is not None:
            path = goal_node.path()
            print(self.goal_path(path))
            return self.goal_path(path)

    def goal_path(self, path):
        """
        Appends the action for each node visited in a solved problem
        @param path: A list of nodes leading to goal state

        @return
        A list of actions performed to reach each new node
        """
        actions = []
        for node in path:
            if node.action is not None:
                actions.append(node.action)
        return actions

    def h(self, node):
        """
        Determines a heuristic value based on the current
        state
        @param node: The current node in the search tree

        @return
        A heuristic value calculation based on a summation of
         manhattan distances for each box to it's closest target space
        """
        state = list(node.state)
        m_dist = 0
        distances = []
        state = state[1]
        for i in range(len(state)-1):
            if state[i] not in self.warehouse.targets:
                for j in range(len(self.warehouse.targets) - 1):
                    box_x, box_y = state[i]
                    target_x, target_y = self.warehouse.targets[j]
                    m_dist = abs(target_x - box_x) + abs(target_y - box_y)
                    distances.append(m_dist)
                m_dist += min(distances)
        return m_dist

    def path_cost(self, c, state1, action, state2):
        """Return the cost of a solution path that arrives at state2 from
        state1 via action, assuming cost c to get up to state1. If the problem
        is such that the path doesn't matter, this function will only look at
        state2.  If the path does matter, it will consider c and maybe state1
        and action. The default method costs 1 for every step in the path."""
        return c + 1

    def result(self, state, action):
        """
        Updates the state of a box and the worker based on the list
        of possible actions.
        @param state: The current state of the warehouse (worker and box locations)
        @param action: List of actions returned by the function actions

        @return
        An updated warehouse state after box and player have been moved
        """
        worker = state[0]
        boxes = state[1]
        box_list = list(boxes)
        Actions = [action]
        for direction in Actions:
            if direction == "Up":
                worker = (worker[0], worker[1]-1)
            elif direction == "Down":
                worker = (worker[0], worker[1]+1)
            elif direction == "Left":
                worker = (worker[0]-1, worker[1])
            elif direction == "Right":
                worker = (worker[0]+1, worker[1])
        #  Player should move and if he moves into a box, the box then needs to move a square in the same direction
            for box in boxes:
                b = box  # temporary variable to hold box coordinates
                if b == worker:
                    if direction == "Up":
                        b = (b[0] , b[1]-1)
                    elif direction == "Down":
                        b = (b[0], b[1]+1)
                    elif direction == "Left":
                        b = (b[0]-1, b[1])
                    elif direction == "Right":
                        b = (b[0]+1, b[1])
                    box_list.remove(box)  # remove the current location of the box
                    box_list.append(b)   # add the location of the new box position
            continue  #  break the loop as the box and/or player has been moved
        return worker,tuple(box_list)  #  updated state

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
    state = (warehouse.worker, tuple(warehouse.boxes))
    puzzle = SokobanPuzzle(warehouse)
    if puzzle.goal_test(state):
        return []
    #sol_ts = search.breadth_first_tree_search(puzzle)
    sol_ts = search.breadth_first_graph_search(puzzle)
    #sol_ts = search.iterative_deepening_search(puzzle)
    #sol_ts = search.depth_first_graph_search(puzzle)
    #sol_ts = search.astar_graph_search(puzzle)
    if sol_ts == None:
        return ['Impossible']
    else:
        ans = puzzle.print_solution(sol_ts)
        return ans

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def can_go_there(warehouse, dst):
    '''
    Determine whether the worker can walk to the cell dst=(row,col)
    without pushing any box.

    @param warehouse: a valid Warehouse object

    @return
      True if the worker can walk to cell dst=(row,col) without pushing any box
      False otherwise

    Logic:
    if the worker is on the y axis
    for the range of the movement to the destination in the y axis, if the x axis meets a box
    then return false, likewise for moving in the x direction return false if box is on same y coordinate
    '''
    X, Y = zip(*warehouse.walls)
    x_size, y_size = 1 + max(X), 1 + max(Y)
    x_movement = []
    for x in range(warehouse.worker[0],dst[0]):  # for the range of the worker to the destination in x axis
        x_movement.append(x)
    y_movement = []
    for y in range(warehouse.worker[1], dst[1]):  # for the range of the worker to the destination in y axis
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
    M = []  # List for macro actions
    macro_wh = warehouse.copy()  #  Copy the warehouse in the original layout
    elem_sol = solve_sokoban_elem(warehouse)  # Solve the warehouse, boxes and worker locations change
    if elem_sol == "Impossible":  # Check if solution was not possible
        return elem_sol
    for action in elem_sol:  # For all the actions in the solution
        boxes = macro_wh.boxes  # boxes = original box locations
        worker = macro_wh.worker
        if action == "Up":
            worker = (macro_wh.worker[0],macro_wh.worker[1]-1)  # update workers location
            for i, box in enumerate(boxes):  # loop the boxes
                if worker == box:  # if worker matches the box location
                     boxes[i] = (box[0], box[1]-1)  # update the box location
                     M.append(((worker[1],worker[0]), "Up"))  # add to the list the new worker (row,col) and the direction
        if action == "Down":
             worker = (macro_wh.worker[0],macro_wh.worker[1]+1)
             for i, box in enumerate(boxes):
                if worker == box:
                    boxes[i] = (box[0], box[1]+1)
                    M.append(((worker[1],worker[0]), "Down"))
        if action == "Left":
            worker = (macro_wh.worker[0]-1,macro_wh.worker[1])
            for i,box in enumerate(boxes):
                if worker == box:
                    boxes[i] = (box[0]-1, box[1])
                    M.append(((worker[1],worker[0]), "Left"))
        if action == "Right":
             worker = (macro_wh.worker[0]+1,macro_wh.worker[1])
             for i,box in enumerate(boxes):
                if worker == box:
                     boxes[i] = (box[0]+1, box[1])
                     M.append(((worker[1],worker[0]), "Right"))
        # for each action save the updated state of the worker and boxes to the copied warehouse
        macro_wh = macro_wh.copy(worker, boxes)
    return M
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def taboo_walls(taboo_tiles, vis, colx):
    """
    This function uses the itertools to find two matching y or x coordinates within
    the list of corner taboos. If a match is found the appropriate helper function is called
    for a match in the y axis check_top_bot_walls() is called where the matching y values will be evaluated
    and determined if the wall is a taboo wall. This logic applies to matching x values and the check_side_walls()
    helper function.

    @param taboo_tiles: A list of all corner taboos
    @param vis: two dimensional array of the warehouse object
    @param colx: the length of the warehouse in the x axis
    """
    def check_top_bot_walls():
            wall_t = 0
            wall_b = 0
            for x in range(x1+1, x2):  # for the x coordinates accompanying the matching y coordinates
                if vis[y1-1][x] == "#" and vis[y1][x] != ".":  # if there is wall and no targets
                    wall_t += 1  #  increase wall counter
                if vis[y1+1][x] == "#" and vis[y1][x] != ".":
                    wall_b += 1
                if wall_t == abs(x2 - (x1+1)):  # if the wall counter matches the difference in x coordinates
                        for wall in range((x1+1), x2):
                            vis[y1][wall] = "X"  # then wall is taboo, update vis array
                if wall_b == abs(x2 - (x1+1)):
                        for w in range((x1+1), x2):
                            vis[y1][w] = "X"
    def check_side_walls():
        wall_l = 0
        wall_r = 0
        for y in range(y1+1, y2):  # for the y coordinates accompanying the matching x coordinates
            if vis[y][x1 - 1] == "#" and vis[y][x1] != ".":  # if there is a wall and no targets
                wall_l += 1  # increase wall counter
            if wall_l == abs(y2 - (y1+1)):  # if wall counter matches the difference in the y direction
                for w in range((y1 + 1), y2):
                   vis[w][x1] = "X"  # update vis
            if x1 + 1 < colx:
                if vis[y][x1 + 1] == "#" and vis[y][x1] != ".":
                    wall_r += 1
                if wall_r == abs(y2 - (y1+1)):
                    for w in range((y1+1), y2):
                        vis[w][x1] = "X"
    for t1, t2 in itertools.permutations(taboo_tiles, 2):
        x1,y1 = t1  # x and y for taboo corner one
        x2,y2 = t2  # x and y for taboo corner two
        if y1 == y2:  # if y1 and y2 match then call helper function
            check_top_bot_walls()
        if x1 == x2:  # if x1 and x2 match then call helper function
            check_side_walls()

def corner_taboos(freetile, warehouse, taboo_tiles, vis):
    """
    Finds and marks all the corners within the warehouse with taboo cells
    @param freetile: List of free tiles on the warehouse object
    @param warehouse: A warehouse object
    @param taboo_tiles: List to hold the taboo tiles found
    @param vis: Two dimensional array of the warehouse
    """
    for tiles in freetile: #Check every single free tile against every possible angle a wall could be
        tile_x, tile_y = tiles  # assign coords to all the blank tiles
        up_x, up_y = tile_x, tile_y - 1  # set up up facing position
        down_x, down_y = tile_x, tile_y + 1  # down ..
        right_x, right_y = tile_x + 1, tile_y  # right..
        left_x, left_y = tile_x - 1, tile_y  # left..
        #  if the worker is facing a corner then append appropriate tile
        if ((up_x, up_y) in warehouse.walls and (right_x, right_y) in warehouse.walls) and vis[tile_y][tile_x] != "." \
                or ((down_x, down_y) in warehouse.walls and (right_x, right_y) in warehouse.walls) \
                and vis[tile_y][tile_x] != "." or ((down_x, down_y) in warehouse.walls and (left_x, left_y) in\
                warehouse.walls) and vis[tile_y][tile_x] != "." or ((up_x, up_y) in warehouse.walls and
                (left_x, left_y) in warehouse.walls) and vis[tile_y][tile_x] != ".":
            taboo_tiles.append(tiles)


def clear_outter_taboos(rowy, colx, vis, warehouse):
    """
    Attempts to clean up taboo cells that are marked outside the warehouse bounds,
    doesn't work perfectly on all warehouse maps.

    @param warehouse: A warehouse object
    @param rowy: Number of rows in the warehouse object
    @param colx: Number of columns in the warehouse object
    @param vis: Two dimensional array of the warehouse
    """
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

def taboo_coordinates(warehouse):
    """
    This function is basically the same as the taboo_cells except it returns the coordinates
    of all the taboo cells instead of a multi-line string

    @param warehouse: A warehouse object

    @:return
    A list of taboo cell coordinates
    """
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if __name__ == "__main__":
    print("main")
    wh = Warehouse()
    wh.read_warehouse_file("./warehouses/warehouse_57.txt")
    print(taboo_cells(wh))
    puz = wh.copy(wh.worker,wh.boxes)
    string = puz.__str__()
    wh.extract_locations(string.split(sep='\n'))
    start = time.time()
    ans = solve_sokoban_elem(wh)
    end = time.time()
    print("Execution Time: ", (end-start))
    print(ans)




