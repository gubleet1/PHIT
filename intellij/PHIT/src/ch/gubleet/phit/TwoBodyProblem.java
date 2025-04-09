package ch.gubleet.phit;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ScheduledFuture;
import java.util.concurrent.TimeUnit;

public class TwoBodyProblem {

    // used for simulation
    private static final double GRAVITATIONAL_CONSTANT = 6.6743e-11;
    private static final double EARTH_WEIGHT = 5.972e24;
    private static final double MOON_WEIGHT = 7.349e22;
    private static final double G_M1 = GRAVITATIONAL_CONSTANT * EARTH_WEIGHT;
    private static final double G_M2 = GRAVITATIONAL_CONSTANT * MOON_WEIGHT;

    // algorithms
    private static final int RUNGE_KUTTA_4 = 0;
    private static final int EULER = 1;

    //initial values
    private static final double earth_x_init = 0;
    private static final double earth_y_init = 0;
    private static final double earth_vx_init = 0;
    private static final double earth_vy_init = 1.131e1;
    private static final double moon_x_init = 3.844e8;
    private static final double moon_y_init = 0;
    private static final double moon_vx_init = 0;
    private static final double moon_vy_init = -9.189e2;

    // --------------------------------------------------------------
    // CONFIGURE VALUES TO ADJUST SIMULATION BEHAVIOR
    // --------------------------------------------------------------
    // 8.64e4 = 1 Day (sufficient for runge kutta 4th order)
    // 3.6e3 = 1 Hour (almost sufficient for euler)
    private static final double TIME_STEP = 8.64e3;
    // RUNGE_KUTTA_4 or EULER
    private static final int ALGORITHM = RUNGE_KUTTA_4;
    // alpha for law of gravitation 1/r^alpha
    private static final double ALPHA = 2.0;

    // only used to calculate reasonable paint refresh rate
    private static final double AVG_MOON_EARTH_DISTANCE = 3.844e8;
    private static final double AVG_MOON_SPEED = 9.189e2;
    private static final int SECONDS_PER_REVOLUTION = 3;
    private static final long STEP_DELAY_MICROSECONDS = Math.round(SECONDS_PER_REVOLUTION * 1e6 / (2.628e6 / TIME_STEP));

    private final Object lock = new Object();
    private ScheduledFuture future;
    private boolean running;

    private List<Position> earthPath;
    private List<Position> moonPath;
    private BufferedImage buffer;

    // earth state
    private double earth_x;
    private double earth_y;
    private double earth_vx;
    private double earth_vy;

    private double earth_dx;
    private double earth_dy;
    private double earth_dvx;
    private double earth_dvy;

    // moon state
    private double moon_x;
    private double moon_y;
    private double moon_vx;
    private double moon_vy;

    private double moon_dx;
    private double moon_dy;
    private double moon_dvx;
    private double moon_dvy;

    private double scale;
    private double origin_x;
    private double origin_y;

    private int paint_d;
    private int unpainted;

    private JPanel cp;
    private SimulationPanel sim;
    private JPanel glassPanel;

    private Color c_earth = new Color(192, 57, 43);
    private Color c_moon = new Color(41, 128, 185);

    private TwoBodyProblem() {
        Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
        double scale_max = Math.max(screen.width, screen.height) / (AVG_MOON_EARTH_DISTANCE * 2);
        paint_d = (int) Math.ceil(10 / (AVG_MOON_SPEED * TIME_STEP * scale_max));
        initialize();
        JFrame frame = new JFrame();
        frame.setTitle("Two Body Problem Simulation");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setResizable(true);
        frame.setMinimumSize(new Dimension(700, 500));
        frame.setExtendedState(JFrame.MAXIMIZED_BOTH);
        createUI();
        frame.setContentPane(cp);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        new TwoBodyProblem();
    }

    private void initialize() {
        earthPath = new LinkedList<>();
        moonPath = new LinkedList<>();
        // earth
        earth_x = earth_x_init;
        earth_y = earth_y_init;
        earth_vx = earth_vx_init;
        earth_vy = earth_vy_init;
        earthPath.add(new Position(earth_x, earth_y));
        // moon
        moon_x = moon_x_init;
        moon_y = moon_y_init;
        moon_vx = moon_vx_init;
        moon_vy = moon_vy_init;
        moonPath.add(new Position(moon_x, moon_y));
        // paint counter
        unpainted = 0;
    }

    @SuppressWarnings("ConstantConditions")
    private void step() {
        // method to solve the differential equations
        if (ALGORITHM == RUNGE_KUTTA_4) {
            rungeKutta4();
        }
        if (ALGORITHM == EULER) {
            euler();
        }
        unpainted++;
        if (unpainted == paint_d) {
            earthPath.add(new Position(earth_x, earth_y));
            moonPath.add(new Position(moon_x, moon_y));
            synchronized (lock) {
                if (buffer != null) {
                    drawLast();
                }
            }
            unpainted = 0;
        }
    }

    private void rungeKutta4() {
        double t = TIME_STEP;
        double[] dt = new double[]{t / 2.0, t / 2.0, t, 0.0};
        double[] c = new double[]{t / 6.0, t / 3.0, t / 3.0, t / 6.0};
        double[] u0 = new double[]{earth_x, earth_y, earth_vx, earth_vy, moon_x, moon_y, moon_vx, moon_vy};
        double[] ut = new double[]{0, 0, 0, 0, 0, 0, 0, 0};
        double[] u = varToU();
        for (int j = 0; j < 4; j++) {
            uToVar(u);
            derivative();
            double[] du = new double[]{earth_dx, earth_dy, earth_dvx, earth_dvy, moon_dx, moon_dy, moon_dvx, moon_dvy};
            for (int i = 0; i < 4; i++) {
                u[i] = u0[i] + dt[j] * du[i];
                ut[i] = ut[i] + c[j] * du[i];
                u[i + 4] = u0[i + 4] + dt[j] * du[i + 4];
                ut[i + 4] = ut[i + 4] + c[j] * du[i + 4];
            }
        }
        for (int i = 0; i < 8; i++) {
            u[i] = u0[i] + ut[i];
        }
        uToVar(u);
    }

    private void euler() {
        double t = TIME_STEP;
        derivative();
        // earth
        earth_x += earth_dx * t;
        earth_y += earth_dy * t;
        earth_vx += earth_dvx * t;
        earth_vy += earth_dvy * t;
        // moon
        moon_x += moon_dx * t;
        moon_y += moon_dy * t;
        moon_vx += moon_dvx * t;
        moon_vy += moon_dvy * t;
    }

    private void derivative() {
        // vector from earth to moon
        double dx = moon_x - earth_x;
        double dy = moon_y - earth_y;
        double l = Math.sqrt(dx * dx + dy * dy);
        // earth
        earth_dx = earth_vx;
        earth_dy = earth_vy;
        double p = ALPHA + 1;
        earth_dvx = (G_M2 / Math.pow(l, p)) * dx;
        earth_dvy = (G_M2 / Math.pow(l, p)) * dy;
        // moon
        moon_dx = moon_vx;
        moon_dy = moon_vy;
        moon_dvx = (-G_M1 / Math.pow(l, p)) * dx;
        moon_dvy = (-G_M1 / Math.pow(l, p)) * dy;
    }

    private void uToVar(double[] u) {
        earth_x = u[0];
        earth_y = u[1];
        earth_vx = u[2];
        earth_vy = u[3];
        moon_x = u[4];
        moon_y = u[5];
        moon_vx = u[6];
        moon_vy = u[7];
    }

    private double[] varToU() {
        double[] u = new double[8];
        u[0] = earth_x;
        u[1] = earth_y;
        u[2] = earth_vx;
        u[3] = earth_vy;
        u[4] = moon_x;
        u[5] = moon_y;
        u[6] = moon_vx;
        u[7] = moon_vy;
        return u;
    }

    private Point2D posToPoint(Position position) {
        double x = origin_x + position.x * scale;
        double y = origin_y + position.y * scale * -1;
        return new Point2D.Double(x, y);
    }

    private void start() {
        if (!running) {
            ScheduledExecutorService scheduler = Executors.newSingleThreadScheduledExecutor();
            future = scheduler.scheduleAtFixedRate(this::step,
                    0, STEP_DELAY_MICROSECONDS, TimeUnit.MICROSECONDS);
            running = true;
        }
    }

    private void stop() {
        if (running) {
            future.cancel(false);
            running = false;
        }
    }

    private void reset() {
        if (running) {
            stop();
        }
        initialize();
        synchronized (lock) {
            buffer = null;
        }
        glassPanel.repaint();
        sim.repaint();
    }

    private void drawLast() {
        drawLast(earthPath);
        drawLast(moonPath);
        glassPanel.repaint();
        sim.repaint();
    }

    private void drawLast(List<Position> path) {
        Graphics2D g = buffer.createGraphics();
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        if (path == earthPath) {
            g.setColor(c_earth);
        }
        if (path == moonPath) {
            g.setColor(c_moon);
        }
        g.setStroke(new BasicStroke(2));
        Point2D from = posToPoint(path.get(path.size() - 1));
        Point2D to = posToPoint(path.get(path.size() - 2));
        g.draw(new Line2D.Double(from, to));
        g.dispose();
    }

    private void drawAll() {
        drawAll(earthPath);
        drawAll(moonPath);
        glassPanel.repaint();
        sim.repaint();
    }

    private void drawAll(List<Position> path) {
        Graphics2D g = buffer.createGraphics();
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        if (path == earthPath) {
            g.setColor(c_earth);
        }
        if (path == moonPath) {
            g.setColor(c_moon);
        }
        g.setStroke(new BasicStroke(2));
        Point2D prev = null;
        for (int i = 0; i < path.size(); i++) {
            Position position = path.get(i);
            Point2D curr = posToPoint(position);
            if (prev == null) {
                prev = curr;
                continue;
            }
            g.draw(new Line2D.Double(prev, curr));
            prev = curr;
        }
        g.dispose();
    }

    private void createUI() {
        cp = new JPanel();
        cp.setLayout(new BorderLayout());
        // simulation
        sim = new SimulationPanel();
        // glass panel
        glassPanel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2d = (Graphics2D) g;
                g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                // earth indicator
                Point2D p = posToPoint(earthPath.get(earthPath.size() - 1));
                Ellipse2D e = new Ellipse2D.Double(p.getX() - 4, p.getY() - 4, 8, 8);
                g.setColor(c_earth);
                g2d.fill(e);
                // moon indicator
                p = posToPoint(moonPath.get(moonPath.size() - 1));
                e = new Ellipse2D.Double(p.getX() - 4, p.getY() - 4, 8, 8);
                g.setColor(c_moon);
                g2d.fill(e);
            }
        };
        glassPanel.setOpaque(false);
        // buttons
        JPanel buttons = new JPanel();
        JButton start = new JButton("Start");
        start.addActionListener(e -> start());
        start.setFocusable(false);
        JButton stop = new JButton("Stop");
        stop.addActionListener(e -> stop());
        stop.setFocusable(false);
        JButton reset = new JButton("Reset");
        reset.addActionListener(e -> reset());
        reset.setFocusable(false);
        buttons.add(start);
        buttons.add(stop);
        buttons.add(reset);
        // content pane
        JLayeredPane layers = new Layers();
        layers.add(glassPanel, Integer.valueOf(1));
        layers.add(sim, Integer.valueOf(0));
        cp.add(layers, BorderLayout.CENTER);
        cp.add(buttons, BorderLayout.SOUTH);
    }

    private static class Position {

        double x;
        double y;

        Position(double x, double y) {
            this.x = x;
            this.y = y;
        }
    }

    private static class Layers extends JLayeredPane {

        Layers() {
            setLayout(new LayoutManager() {
                @Override
                public void addLayoutComponent(String name, Component comp) {
                    // nothing
                }

                @Override
                public void removeLayoutComponent(Component comp) {
                    // nothing
                }

                @Override
                public Dimension preferredLayoutSize(Container parent) {
                    return new Dimension(Integer.MAX_VALUE, Integer.MAX_VALUE);
                }

                @Override
                public Dimension minimumLayoutSize(Container parent) {
                    return new Dimension(0, 0);
                }

                @Override
                public void layoutContainer(Container parent) {
                    for (int i = 0; i < parent.getComponentCount(); i++) {
                        parent.getComponent(i).setBounds(parent.getBounds());
                    }
                }
            });
        }
    }

    private class SimulationPanel extends JPanel {

        SimulationPanel() {
            setBorder(BorderFactory.createLineBorder(Color.BLACK));
            setBackground(Color.WHITE);
        }

        private BufferedImage getBuffer() {
            if (buffer == null) {
                Insets i = getInsets();
                int w = Math.max(1, getWidth() - i.left - i.right);
                int h = Math.max(1, getHeight() - i.top - i.bottom);
                buffer = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
                origin_x = i.left + (w / 2.0);
                origin_y = i.top + (h / 2.0);
                double scale_w = w / (AVG_MOON_EARTH_DISTANCE * 2);
                double scale_h = h / (AVG_MOON_EARTH_DISTANCE * 2);
                scale = Math.min(scale_w, scale_h);
                drawAll();
            }
            return buffer;
        }

        @Override
        public void invalidate() {
            synchronized (lock) {
                buffer = null;
                super.invalidate();
            }
        }

        @Override
        protected void paintComponent(Graphics g) {
            synchronized (lock) {
                super.paintComponent(g);
                Insets i = getInsets();
                g.drawImage(getBuffer(), i.left, i.top, this);
            }
        }
    }
}
